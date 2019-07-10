/* calc_signal.c -- based on m3d2s.f by I-Yang Lee
 * Karin Lagergren
 *
 * This module contains the main interface to the signal calculation
 * code. 
 *
 * To use: 
 * -- call signal_calc_init. This will initialize geometry, fields,
 *       drift velocities etc.
 * -- call get_signal
 *
 */
/* TODO: see FIXME's below
   charge_trapping is just a placeholder ATM. Should it be defined
   in the fields module?
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "mjd_siggen.h"
#include "calc_signal.h"
#include "point.h"
#include "detector_geometry.h"
#include "fields.h"

#define HOLE_CHARGE 1.0
#define ELECTRON_CHARGE -1.0
/* the following is the diffusion coefficient for holes in Ge at 77K
               at low field (~ 100 V/cm) 
   the diffusion coefficient drops at higher fields, and higher temperatures
   see Jacoboni et al., Phys. Rev. B24, 2 (1981) 1014-1026.
   size sigma = sqrt(2Dt), t = time, D = mu*k*T/e
   mu = mobility, k = Boltzmann const., T = temp, e = electron charge
   mu_h = 4e4 cm^2/V/s, mu_e = 5e4 cm^2/V/s at 77K, so
   D_h = 265 cm^2/s, D_e = 332 cm^2/s
   and goes down roughly as 1/Temp (since mu goes as T^-1.7 or T^-2.3)

   we also convert (2Dt) from sigma-squared to FWHM-squared

   for Si at 300K, 
   mu_h = 450 cm^2/V/s, mu_e = 1500 cm^2/V/s, so
   D_h = 12 cm^2/s, D_e = 39 cm^2/s
*/
/*  here are some definitions used for an old method, where I calculated FWHM-squared:
// germanium:  2Dt = sigma^2; we want  FWHM^2 in mm^2 / ns
#define TWO_TIMES_DIFFUSION_COEF_H \
        (2.0 * 2.355*2.355 * 2.65e-5 * setup->step_time_calc * 77.0/setup->xtal_temp)
#define TWO_TIMES_DIFFUSION_COEF_E \
        (2.0 * 2.355*2.355 * 3.32e-5 * setup->step_time_calc * 77.0/setup->xtal_temp)
// silicon:
#define TWO_TIMES_DIFFUSION_COEF_H_Si \
        (1.3e-5 * setup->step_time_calc * 300.0/setup->xtal_temp)
#define TWO_TIMES_DIFFUSION_COEF_E_Si \
        (4.3e-5 * setup->step_time_calc * 300.0/setup->xtal_temp)
*/
/* In the new method, I use dsigma/dt = D/sigma to calculate FWHM */
#define DIFFUSION_COEF   (setup->v_over_E * 0.67)
/* above is my own approximate parameterization of measurements of Jacoboni et al.
   0.67 = 2.355 * 2.355 * 0.12    to get D in mm2/s, and scaled to FWHM2/sigma2
   v_over_E = drift velocity / electric field   ~  mu
   note that Einstein's equation is D = mu*kT/e
   kT/e ~ 0.007/V ~ 0.07 mm/Vcm, => close enough to 0.12, okay
 */

/* prototypes for module-private functions*/
//static int make_signal(point pt, float *signal, float q, MJD_Siggen_Setup *setup);
//static float charge_trapping(vector dx, float q); //trapping

/* signal_calc_init
   read setup from configuration file,
   then read the electric field and weighting potential,
   and initialize the signal calculation variables
   returns 0 for success
*/
int signal_calc_init(char *config_file_name, MJD_Siggen_Setup *setup) {

  if (read_config(config_file_name, setup)) return 1;

  TELL_CHATTY("r: %.2f  z: %.2f\n", setup->xtal_radius, setup->xtal_length);
  setup->ntsteps_out = setup->time_steps_calc /
    lrintf(setup->step_time_out/setup->step_time_calc);
  TELL_NORMAL("Will use %d time steps in calculations, each %.2f ns long;\n"
	      "the output signals will have %d time steps, each %.2f ns long\n", 
	      setup->time_steps_calc, setup->step_time_calc,
	      setup->ntsteps_out, setup->step_time_out);

  TELL_NORMAL("Reading field data...\n");
  if (field_setup(setup) != 0) return -1;
  
  if ((setup->dpath_e = (point *) malloc(setup->time_steps_calc*sizeof(point))) == NULL ||
      (setup->dpath_h = (point *) malloc(setup->time_steps_calc*sizeof(point))) == NULL) {
    error("Path malloc failed\n");
    return -1;
  }

  tell("Setup of signal calculation done\n");
  return 0;
}

/* get_signal
   calculate signal for point pt. Result is placed in signal_out array
   returns -1 if outside crystal
   if signal_out == NULL => no signal is stored
*/
int get_signal(point pt, float *signal_out, MJD_Siggen_Setup *setup) {
  static float *signal, *sum, *tmp;
  static int tsteps = 0;
  float w, x, y;
  char  tmpstr[MAX_LINE];
  int   j, k, l, dt, err, comp_f;

  /* first time -- allocate signal and sum arrays */
  if (tsteps != setup->time_steps_calc) {
    tsteps = setup->time_steps_calc;
    if ((signal = (float *) malloc(tsteps*sizeof(*signal))) == NULL ||
	(tmp    = (float *) malloc(tsteps*sizeof(*tmp))) == NULL ||
	(sum    = (float *) malloc(tsteps*sizeof(*sum))) == NULL) {
      error("malloc failed in get_signal\n");
      return -1;
    }
  }

  for (j = 0; j < tsteps; j++) signal[j] = 0.0;

  if (outside_detector(pt, setup)) {
    TELL_CHATTY("Point %s is outside detector!\n", pt_to_str(tmpstr, MAX_LINE, pt));
    return -1;
  }
  TELL_CHATTY("Calculating signal for %s...\n", pt_to_str(tmpstr, MAX_LINE, pt));

  memset(setup->dpath_e, 0, tsteps*sizeof(point));
  memset(setup->dpath_h, 0, tsteps*sizeof(point));

  err = make_signal(pt, signal, ELECTRON_CHARGE, setup);
  err = make_signal(pt, signal, HOLE_CHARGE, setup);
  /* make_signal returns 0 for success; require hole signal but not electron */

  /* change from current signal to charge signal, i.e.
     each time step contains the summed signals of all previous time steps */
  for (j = 1; j < tsteps; j++) signal[j] += signal[j-1];

  if (signal_out != NULL) {

    if (setup->charge_cloud_size > 0.001 || setup->use_diffusion) {
      /* convolute with a Gaussian to correct for charge cloud size
	 and initial velocity
	 charge_cloud_size = initial FWHM of charge cloud, in mm,
	 NOTE this uses initial velocity of holes only;
	 this may not be quite right if electron signal is strong */
      /* difference in time between center and edge of charge cloud */
      dt = (int) (1.5f + setup->charge_cloud_size /
		  (setup->step_time_calc * setup->initial_vel));
      if (setup->initial_vel < 0.00001f) dt = 0;
      TELL_CHATTY("Initial vel, size, dt = %f mm/ns, %f mm, %d steps\n",
		  setup->initial_vel, setup->charge_cloud_size, dt);
      if (setup->use_diffusion) {
	dt = (int) (1.5f + setup->final_charge_size /
		    (setup->step_time_calc * setup->final_vel));
	TELL_CHATTY("  Final vel, size, dt = %f mm/ns, %f mm, %d steps\n",
		    setup->final_vel, setup->final_charge_size, dt);
      }
      if (dt > 1) {
	/* Gaussian */
	w = ((float) dt) / 2.355;
	l = dt/10;     // use l to speed up convolution of waveform with gaussian;
	if (l < 1) {   // instead of using every 1-ns step, use steps of FWHM/10
	  l = 1;
	} else if (setup->step_time_out > setup->preamp_tau) {
	  if (l > setup->step_time_out/setup->step_time_calc)
	    l = setup->step_time_out/setup->step_time_calc;
	} else {
	  if (l > setup->preamp_tau/setup->step_time_calc)
	    l = setup->preamp_tau/setup->step_time_calc;
	}
	// TELL_CHATTY(">> l: %d\n", l);
	for (j = 0; j < tsteps; j++) {
	  sum[j] = 1.0;
	  tmp[j] = signal[j];
	}
	for (k = l; k < 2*dt; k+=l) {
	  x = ((float) k)/w;
	  y = exp(-x*x/2.0);
	  for (j = 0; j < tsteps - k; j++){
	    sum[j] += y;
	    tmp[j] += signal[j+k] * y;
	    sum[j+k] += y;
	    tmp[j+k] += signal[j] * y;
	  }
	}
        for (j = 0; j < tsteps; j++){
          signal[j] = tmp[j]/sum[j];
        }
      }
    }

    /* now, compress the signal and place it in the signal_out array;
       truncate the signal if time_steps_calc % ntsteps_out != 0 */
    comp_f = setup->time_steps_calc/setup->ntsteps_out;
    for (j = 0; j < setup->ntsteps_out; j++) signal_out[j] = 0;
    for (j = 0; j < setup->ntsteps_out*comp_f; j++)
      signal_out[j/comp_f] += signal[j]/comp_f;

    /* do RC integration for preamp risetime */
    if (setup->preamp_tau/setup->step_time_out >= 0.1f)
      rc_integrate(signal_out, signal_out,
		   setup->preamp_tau/setup->step_time_out, setup->ntsteps_out);
  }

  /* make_signal returns 0 for success; require hole signal but not electron */
  if (err) return -1;
  return 1;
}

/* make_signal
   Generates the signal originating at point pt, for charge q
   returns 0 for success
*/
int make_signal(point pt, float *signal, float q, MJD_Siggen_Setup *setup) {
  static float wpot, dwpot;
  char   tmpstr[MAX_LINE];
  point  new_pt;
  vector v, dx;
  float  vel0, vel1 = 0, wpot_old=-1;
  // double diffusion_coeff;
  double repulsion_fact = 0.0, ds2, ds3, dv, ds_dt;
  int    ntsteps, i, t, n, collect2pc, low_field=0, surface_drift=0;

  new_pt = pt;
  collect2pc = ((q > 0 && setup->impurity_z0 < 0) ||  // holes for p-type 
		(q < 0 && setup->impurity_z0 > 0));   // electrons for n-type
  /*
  if (q > 0) {
    diffusion_coeff = TWO_TIMES_DIFFUSION_COEF_H;
  } else {
    diffusion_coeff = TWO_TIMES_DIFFUSION_COEF_E;
  }
  */
  ntsteps = setup->time_steps_calc;
  for (t = 0; drift_velocity(new_pt, q, &v, setup) >= 0; t++) { 
    if (q > 0) {
      setup->dpath_h[t] = new_pt;
    } else {
      setup->dpath_e[t] = new_pt;
    }
    if (collect2pc) {
      if (t == 0) {
	vel1 = setup->final_vel = setup->initial_vel = vector_length(v);
	setup->final_charge_size = setup->charge_cloud_size;
	if (setup->use_diffusion) {
	  if (setup->final_charge_size < 0.01) setup->final_charge_size = 0.01;
	  /* for a spherically symmetric charge cloud, the equivalent
	     delta-E at a distance of 1 sigma from the cloud center is
	     dE = Q/(4*pi*epsilon*sigma^2)  (Q is charge inside the 3D 1-sigma envelope)
	     dE (V/cm) = Q (C) * 1/(4*pi*epsilon) (N m2 / C2) / sigma2 (mm2)
	     1 V/m = 1 N/C
	     dE (V/cm) = Q (C) * 1/(4*pi*epsilon) (V m / C) / sigma2 (mm2)
	     dE (V/cm) = repulsion_fact * FWHM/sigma / (FWHM^2) (mm2), so
	     repulsion_fact = (FWHM/sigma)^3 * Q (C) * 1/(4*pi*epsilon) (V m / C) * mm/m * mm/cm
	  */
	  if (setup->energy > 0.1) {  // set up charge cloud self-repulsion
	    repulsion_fact = setup->energy * 0.67*0.67*0.67 / 0.003; // charge in 1 sigma (3D)
	    repulsion_fact /= 6.241e18;        // convert to Coulombs
	    repulsion_fact *= 9.0e13/16.0;     // 1/(4*pi*epsilon)  (N m2 / C2) * 1e4
	    repulsion_fact *= 2.355*2.355*2.355;      // convert FWHM to sigma
	  }
	}
	TELL_CHATTY("initial v: %f (%e %e %e)\n",
		    setup->initial_vel, v.x, v.y, v.z);
      } else if (setup->use_diffusion) {
	vel0 = vel1;
	vel1 = vector_length(v);
	setup->final_charge_size *= vel1/vel0;  // effect of acceleration
	// include effects of acceleration and diffusion on cloud size
	dv = repulsion_fact * setup->dv_dE /        // effect of repulsion
	        (setup->final_charge_size*setup->final_charge_size);
	// FIXME? this next line could more more fine-grained
	if (dv > 0.05) dv = 0.05;  // on account of drift velocity saturation
	ds_dt = dv + DIFFUSION_COEF/setup->final_charge_size;  // effect of diffusion
	if (ds_dt > 0.05 || ds_dt * setup->step_time_calc > 0.1) {
	  // nonlinear growth due to small size; need more careful calculation
	  TELL_CHATTY("ds_dt = %.2f; size = %.2f", ds_dt, setup->final_charge_size);
	  // ds_dt = 0.05;  // artificially limit nonlinear growth
	  ds2 = 2.0 * DIFFUSION_COEF * setup->step_time_calc; // increase^2 from diff.
	  ds3 = (setup->final_charge_size*setup->final_charge_size *
		 (setup->final_charge_size +
		  3.0 * dv * setup->step_time_calc));         // FWHM^3 after repulsion
	  setup->final_charge_size = sqrt(ds2 + pow(ds3, 0.6667)); 
	  TELL_CHATTY(" -> %.2f\n", setup->final_charge_size);
	} else {
	  setup->final_charge_size +=  ds_dt * setup->step_time_calc;  // effect of diff. + rep.
	}
      }
    }

    TELL_CHATTY("pt: (%.2f %.2f %.2f), v: (%e %e %e)",
		new_pt.x, new_pt.y, new_pt.z, v.x, v.y, v.z);
    if (t >= ntsteps - 2) {
      if (collect2pc || wpot > WP_THRESH_ELECTRONS) {
	/* for p-type, this is hole or electron+high wp */
	TELL_CHATTY("\nExceeded maximum number of time steps (%d)\n", ntsteps);
	low_field = 1;
	// return -1;
      }
      break;
    }
    if (wpotential(new_pt, &wpot, setup) != 0) {
      TELL_NORMAL("\nCan calculate velocity but not WP at %s!\n",
		  pt_to_str(tmpstr, MAX_LINE, new_pt));
      return -1;
    }
    if (wpot < 0.0) wpot = 0.0;
    TELL_CHATTY(" -> wp: %.4f\n", wpot);
    if (t > 0) signal[t] += q*(wpot - wpot_old);
    // FIXME? Hack added by DCR to deal with undepleted point contact
    if (wpot >= 0.999 && (wpot - wpot_old) < 0.0002) {
      low_field = 1;
      break;
    }
    wpot_old = wpot;

    dx = vector_scale(v, setup->step_time_calc);
    if (surface_drift && dx.z < 0) {
      dx.x *= setup->surface_drift_vel_factor;  // Hmmm... should the default be zero or one?
      dx.y *= setup->surface_drift_vel_factor;
      dx.z = 0;
    }
    new_pt = vector_add(new_pt, dx);
    // q = charge_trapping(dx, q); //FIXME

    // look for charges on passivated surface of a PPC detector
    if (new_pt.z < 0 &&                                    // at or below surface, and
        (setup->wrap_around_radius <= setup->pc_radius ||  // this is a PPC detector
         new_pt.x*new_pt.x + new_pt.y*new_pt.y <           // or point is inside wrap-around
         setup->wrap_around_radius*setup->wrap_around_radius)) {
      TELL_CHATTY(" -> Passivated surface! q = %.2f  r = %.2f\n",
                  q, sqrt(new_pt.x*new_pt.x + new_pt.y*new_pt.y));
      //break;
      surface_drift = 1;
      new_pt.z = 0;
    }

  }
  if (t == 0) {
    TELL_CHATTY("The starting point %s is outside the field.\n",
		pt_to_str(tmpstr, MAX_LINE, pt));
    return -1;
  }

  if (low_field) {
    TELL_CHATTY("Too many time steps or low field; this may or may not be a problem.\n");
  } else {
    TELL_CHATTY("Drifted to edge of field grid, point: %s q: %.2f\n", 
		pt_to_str(tmpstr, MAX_LINE, new_pt), q);

    /* now we are outside the electric grid. figure out how much we must
       drift to get to the crystal boundary */
    for (n = 0; n+t < ntsteps; n++){
      new_pt = vector_add(new_pt, dx);
      if (q > 0) setup->dpath_h[t+n] = new_pt;
      else setup->dpath_e[t+n] = new_pt;
      if (outside_detector(new_pt, setup)) break;
    }
    if (n == 0) n = 1; /* always drift at least one more step */
    // TELL_CHATTY(
    TELL_NORMAL("q: %.1f t: %d n: %d ((%.2f %.2f %.2f)=>(%.2f %.2f %.2f))\n", 
		q, t, n, pt.x, pt.y, pt.z, new_pt.x, new_pt.y, new_pt.z);

    if (n + t >= ntsteps){
      if (q > 0 || wpot > WP_THRESH_ELECTRONS) { /* hole or electron+high wp */
	TELL_CHATTY("Exceeded maximum number of time steps (%d)\n", ntsteps);
	return -1;  /* FIXME DCR: does this happen? could this be improved? */
      }
      n = ntsteps -t;
    }
    /* make WP go gradually to 1 or 0 */
    if (wpot > 0.3) {
      dwpot = (1.0 - wpot)/n;
    } else {
      dwpot = - wpot/n;
    }

    /*now drift the final n steps*/
    dx = vector_scale(v, setup->step_time_calc);
    if (new_pt.z > 0) {               // charges NOT on passivated surface
      for (i = 0; i < n; i++){
        signal[i+t] += q*dwpot;
        // q = charge_trapping(dx, q); //FIXME
      }
    }
  }
  TELL_CHATTY("q:%.2f pt: %s\n", q, pt_to_str(tmpstr, MAX_LINE, pt));
  if (q > 0) setup->final_vel = vector_length(v);

  return 0;
}

//FIXME -- placeholder function. Even parameter list is dodgy
/*
static float charge_trapping(vector dx, float q){
  return q;
}
*/

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
int signal_calc_finalize(MJD_Siggen_Setup *setup){
  fields_finalize(setup);
  free(setup->dpath_h);
  free(setup->dpath_e);
  return 0;
}

int drift_path_e(point **pp, MJD_Siggen_Setup *setup){
  *pp = setup->dpath_e;
  return setup->time_steps_calc;
}
int drift_path_h(point **pp, MJD_Siggen_Setup *setup){
  *pp = setup->dpath_h;
  return setup->time_steps_calc;
}

/* tell
   write to stdout, provided that verb_level is above the threshold */
void tell(const char *format, ...){
  va_list ap;

  va_start(ap, format);
  vprintf(format, ap);
  va_end(ap);
  return;
}

/*error
  report error messages to stderr */
void error(const char *format, ...) {
  va_list ap;

  va_start(ap, format);
  vfprintf(stderr, format, ap);
  va_end(ap);
  return;
}
