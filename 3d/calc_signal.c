 /* calc_signal.c -- based on m3d2s.f by I-Yang Lee
 * Karin Lagergren
 *
 * This module contains the main interface to the signal calculation
 * code. 
 *
 * To use: 
 * -- call signal_calc_init. This will initialize geometry, fields,
 *       drift velocities etc.
 * -- call hit_segment/get_signal
 *
 *  David Radford   Oct 2019
 *  Updated code by Karin Lagergren to match new config file arrangement and new fieldgen
 *
 */
/* TODO: see FIXME's below
   charge_trapping is just a placeholder ATM. Should it be defined
   in the fields module?
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mjd_siggen.h"
#include "calc_signal.h"
#include "point.h"
#include "detector_geometry.h"
#include "fields.h"

#define HOLE_CHARGE 1.0
#define ELECTRON_CHARGE -1.0

/* prototypes for module-private functions */
static int   make_signal(MJD_Siggen_Setup *setup, point pt, float *signal, float q);
static int   sum_signal(MJD_Siggen_Setup *setup, float *s);         //modifies s
static float charge_trapping(point pt, float distance, float q);
static int   zero_signal(MJD_Siggen_Setup *setup, float *signal);
int   rc_integrate(float *s_in, float *s_out, float tau, int time_steps);


/* signal_calc_init
   read setup from configuration file,
   then read the electric field and weighting potentials,
   and initialize the signal calculation variables
   returns 0 for success
*/
int signal_calc_init(MJD_Siggen_Setup *setup)
{
  char fname[256];

  strncpy(fname, setup->config_file_name, 256);
  tell(TERSE, "Reading configuration file: %s\n", fname);
  if (read_config(fname, setup)) {
    error("Read config failed\n");
    return -1;
  }
  strncpy(setup->config_file_name, fname, 256);
  setup->ntsteps_out = setup->time_steps_calc * setup->step_time_calc / setup->step_time_out;

  tell(TERSE, "Reading fields\n");
  if (field_setup(setup) != 0) {
    error("Field setup failed\n");
    return -1;
  }  
  if ((setup->dpath_e = malloc(setup->time_steps_calc * sizeof(*setup->dpath_e))) == NULL ||
      (setup->dpath_h = malloc(setup->time_steps_calc * sizeof(*setup->dpath_h))) == NULL) {
    error("Malloc failed\n");
    return -1;
  }

  tell(TERSE, "Setup of signal calculation done\n");
  return 0;
}

/* get_signal
   calculate signal for point pt. Result is placed in signal_out array
   which is assumed to have the appropriate size (nsegments * setup->ntsteps_out)
   returns segment number or -1 if outside crystal
   if signal_out == NULL => no signal is stored
*/
int get_signal(MJD_Siggen_Setup *setup, point pt, float *signal_out)
{
  static float *signal;
  static int tsteps = 0;
  char tmpstr[MAX_LINE];
  int  j, comp_f;

  /*first time -- allocate memory for signal array */
  if (tsteps != setup->time_steps_calc) {
    if ((signal = malloc(setup->time_steps_calc * sizeof(*signal))) == NULL) return -1;
    tsteps = setup->time_steps_calc;
  }

  pt_to_str(tmpstr, MAX_LINE, pt);
  if (!in_crystal(setup, pt)) return -1;
  tell(CHATTY, "point %s is in crystal\n", tmpstr);

  zero_signal(setup, signal);
  memset(setup->dpath_e, 0, setup->time_steps_calc*sizeof(*setup->dpath_e));
  memset(setup->dpath_h, 0, setup->time_steps_calc*sizeof(*setup->dpath_h));

  tell(CHATTY, " @@@@@ Signal for %s\n", tmpstr);
  if (make_signal(setup, pt, signal, ELECTRON_CHARGE) ||
      make_signal(setup, pt, signal, HOLE_CHARGE)) return -1;
  /* make_signal returns 0 for success */
  sum_signal(setup, signal);

  if (signal_out != NULL) {
    /* now, compress the signal and place it in the signal_out array */
    comp_f = setup->time_steps_calc/setup->ntsteps_out;
    for (j = 0; j < setup->ntsteps_out; j++) signal_out[j] = 0;
    /* truncate the signal if setup->time_steps_calc % setup->ntsteps_out != 0 */
    for (j = 0; j < setup->ntsteps_out*comp_f; j++)
      signal_out[j/comp_f] += signal[j]/comp_f;
  }

  /* do RC integration for preamp risetime */
  if (setup->preamp_tau/setup->step_time_out >= 0.1f)
    rc_integrate(signal_out, signal_out,
                 setup->preamp_tau/setup->step_time_out, setup->ntsteps_out);
  return 1;
}


static int zero_signal(MJD_Siggen_Setup *setup, float *signal)
{
  int j;
  
  for (j = 0; j < setup->time_steps_calc; j++) signal[j] = 0.0;
  return 0;
}


/* make_signal
   Generates the signal originating at point pt, for charge q
   returns 0 for success
*/
static int make_signal(MJD_Siggen_Setup *setup, point pt, float *signal, float q)
{
  static float wpot, wpot2, wpot_old, dwpot;

  char tmpstr[MAX_LINE];
  point new_pt, prev_pt;
  vector v, dx;
  float dist;
  int t, n, i, vr, stop_drifting = 0;


  prev_pt = new_pt = pt;
  for (t = 0; ((vr = drift_velocity(setup, new_pt, q, &v)) >= 0 && !stop_drifting); t++) { 
    if (q > 0) setup->dpath_h[t] = new_pt;
    else setup->dpath_e[t] = new_pt;
    tell(CHATTY, "t: %d  pt: (%.2f %.2f %.2f)\n", t, new_pt.x,new_pt.y, new_pt.z);
    tell(CHATTY, "v: (%e %e %e)\n", v.x, v.y, v.z);
    if (wpotential(setup, new_pt, &wpot) != 0) {
      pt_to_str(tmpstr, MAX_LINE, new_pt);
      tell(NORMAL,
	   "Can calculate velocity but not weighting potentials at %s!\n",
	   tmpstr);
      return -1;
    }

   if ((q * setup->impurity_z0 < 0 &&    // drifting to point contact  (e in n-type or h in p-type)
         1.0-wpot < 5.0e-5) ||
        (q * setup->impurity_z0 > 0 &&    // drifting to outside (h in n-type or e in p-type)
         wpot< 5.0e-5)) {
      pt_to_str(tmpstr, MAX_LINE, new_pt);
      tell(CHATTY, "Reached full WP at %s; WP = %9.6f\n", tmpstr, wpot);
      stop_drifting = 2;
    }
    /* relax the limit on this check if we have reached a boundary
       i.e field extrapolation is require) */
    if (vr > 0 && new_pt.z < setup->xtal_length - setup->xtal_grid) {  // CHECKME
      if ((q * setup->impurity_z0 < 0 &&      // drifting to point contact (e in n-type or h in p-type)
           1.0-wpot < 0.05 &&
           wpot - wpot_old > 0.0) ||  // make sure local drift is in correct direction
          (q * setup->impurity_z0 > 0 &&      // drifting to outside (h in n-type or e in p-type)
           wpot < 0.05 &&
           wpot - wpot_old < 0.0)) {  // make sure local drift is in correct direction
        pt_to_str(tmpstr, MAX_LINE, new_pt);
        tell(CHATTY, "Reached boundary at %s; WP = %9.6f\n", tmpstr, wpot);
        stop_drifting = 1;
      }
    }
    if (t >= setup->time_steps_calc - 2) stop_drifting = 1; // have run out of time...

    if (t > 0) signal[t] += q*(wpot - wpot_old);
    wpot_old = wpot;
    dx = vector_scale(v, setup->step_time_calc);
    prev_pt = new_pt;
    new_pt = vector_add(new_pt, dx);
    dist = vector_length(dx);
    q = charge_trapping(new_pt, dist, q);   // FIXME
  }

  if (t == 0) {
    pt_to_str(tmpstr, MAX_LINE, pt);
    tell(CHATTY, "The starting point %s is outside the field.\n", tmpstr);
    return -1;
  }
  /*check if we have drifted out of back of detector (not to contact )*/
  if (new_pt.z > setup->xtal_length) {
    tell(CHATTY, "Drifted out of back end of detector.\n");
    // return -1;
    new_pt.z = setup->xtal_length; // FIXME: change horizontal velocity??
  }

  pt_to_str(tmpstr, MAX_LINE, new_pt);
  tell(CHATTY, "Drifted to edge of field grid, point: %s q: %.2f  sd: %d\n", 
       tmpstr, q, stop_drifting);

  /* now we are outside the electric field grid
     decide whether we need to keep drifting to make WPs go to zero */
  if (stop_drifting) return 0;

  /* figure out how much we must drift to get to the crystal boundary */
  for (n = 0; in_crystal(setup, new_pt) &&  n+t < setup->time_steps_calc; n++) {
    new_pt = vector_add(new_pt, dx);
    if (q > 0) setup->dpath_h[t+n] = new_pt;
    else setup->dpath_e[t+n] = new_pt;
    if (n * setup->step_time_calc > 50 &&  // final drift is longer than 50 ns
        new_pt.z < setup->xtal_length)     // not stuck on the passivated surface
      break;                               // no extra steps beyond 50 ns
    
  }
  if (n == 0) n = 1; /* always drift at least one more step */

  tell(CHATTY, "q: %.1f t: %d n: %d ((%.2f %.2f %.2f)=>(%.2f %.2f %.2f))\n", 
       q,t,n, pt.x, pt.y, pt.z, new_pt.x, new_pt.y, new_pt.z);

  if (n + t >= setup->time_steps_calc) {
    tell(CHATTY, "Exceeded maximum number of time steps (%d)\n", setup->time_steps_calc);
    /* check WP's to see if we have produced most of the signal */
    if (wpot < 0.95 && wpotential(setup, new_pt, &wpot2) != 0) {
      tell(CHATTY, "Cannot finish drifting to make at least 95\% of signal.\n");
      return -1;  /* FIXME DCR: could this be improved? */
    }
    /* drift to new_pt and wpot2 */
      dwpot = (wpot2 - wpot)/n;
  } else {
    /* weighting pot. is 1 at edge for hit segment, 0 for other segments.
       Make it so, gradually if applicable */
    dwpot = (1.0 - wpot)/n;
  }
  if (dwpot <  0) {
    tell(CHATTY, "Cannot complete drifting; WP, dWP: %7.4f %7.4f\n",
         wpot, dwpot);
    return -1;  /* FIXME DCR: could this be improved? */
  }

  /* now drift the final n steps */
  tell(CHATTY, " >>> completing drift: t, n = %4d, %4d (%4d);"
       "  WP, dWP: %7.4f %7.4f\n",
       t, n, t+n, wpot, dwpot);
  dx = vector_scale(v, setup->step_time_calc);
  dist = vector_length(dx);
  for (i = 1; i <= n; i++) {
    signal[i+t-1] += q*dwpot;
    q = charge_trapping(prev_pt,dist, q);   // FIXME    
  }

  pt_to_str(tmpstr, MAX_LINE, pt);
  tell(CHATTY, "q:%.2f pt: %s\n", q, tmpstr);

  return 0;
}

/* modifies s. each time step will contain the summed signals of 
   all previous time steps */
static int sum_signal(MJD_Siggen_Setup *setup, float *s)
{
  int j;

  for (j = 1; j < setup->time_steps_calc; j++) {
    s[j] += s[j-1];
  }
  return 0;
}

// FIXME -- placeholder function. Even parameter list is dodgy
static float charge_trapping(point pt, float distance, float q)
{
  return q;
}

int rc_integrate(float *s_in, float *s_out, float tau, int time_steps){
  int   j;
  float s_in_old, s;  /* DCR: added so that it's okay to
			 call this function with s_out == s_in */
  
  if (tau < 1.0f) {
    for (j = time_steps-1; j > 0; j--) s_out[j] = s_in[j-1];
    s_out[0] = 0.0;
  } else {
    s_in_old = s_in[0];
    s_out[0] = 0.0;
    for (j = 1; j < time_steps; j++) {
      s = s_out[j-1] + (s_in_old - s_out[j-1])/tau;
      s_in_old = s_in[j];
      s_out[j] = s;
    }
  }
  return 0;
}

/* signal_calc_finalize
 * Clean up (free arrays, close open files...)
 */
int signal_calc_finalize(MJD_Siggen_Setup *setup)
{
  fields_finalize(setup);
  free(setup->dpath_h);
  free(setup->dpath_e);
  return 0;
}

int drift_path_e(MJD_Siggen_Setup *setup, point **pp)
{
  *pp = setup->dpath_e;
  return setup->time_steps_calc;
}
int drift_path_h(MJD_Siggen_Setup *setup ,point **pp)
{
  *pp = setup->dpath_h;
  return setup->time_steps_calc;
}
