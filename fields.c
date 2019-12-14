/* fields_ppc.c -- based on m3d2s.f by I-Yang Lee 
 * Karin Lagergren
 *
 * This module handles the electric field and weighting potential and 
 * calculates drift velocities
 *
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "mjd_siggen.h"
#include "point.h"
#include "cyl_point.h"
#include "fields.h"
#include "detector_geometry.h"
#include "calc_signal.h"

#define MAX_FNAME_LEN 512

static int nearest_field_grid_index(cyl_pt pt, cyl_int_pt *ipt, MJD_Siggen_Setup *setup);
static int grid_weights(cyl_pt pt, cyl_int_pt ipt, float out[2][2], MJD_Siggen_Setup *setup);
static cyl_pt efield(cyl_pt pt, cyl_int_pt ipt, MJD_Siggen_Setup *setup);
static int setup_efield(MJD_Siggen_Setup *setup);
static int setup_wp(MJD_Siggen_Setup *setup);
static int setup_velo(MJD_Siggen_Setup *setup);
static int efield_exists(cyl_pt pt, MJD_Siggen_Setup *setup);

/* field_setup
   given a field directory file, read electic field and weighting
   potential tables from files listed in directory
   returns 0 for success
*/
int field_setup(MJD_Siggen_Setup *setup){

  setup->rmin  = 0;
  setup->rmax  = setup->xtal_radius;
  setup->rstep = setup->xtal_grid;
  setup->zmin  = 0;
  setup->zmax  = setup->xtal_length;
  setup->zstep = setup->xtal_grid;
  if (setup->xtal_temp < MIN_TEMP) setup->xtal_temp = MIN_TEMP;
  if (setup->xtal_temp > MAX_TEMP) setup->xtal_temp = MAX_TEMP;

  TELL_NORMAL("rmin: %.2f rmax: %.2f, rstep: %.2f\n"
	      "zmin: %.2f zmax: %.2f, zstep: %.2f\n"
	      "Detector temperature is set to %.1f K\n",
	      setup->rmin, setup->rmax, setup->rstep,
	      setup->zmin, setup->zmax, setup->zstep,
	      setup->xtal_temp);

  if (setup_velo(setup) != 0){
    error("Failed to read drift velocity data from file: %s\n", 
	  setup->drift_name);
    return -1;
  }
  if (setup_efield(setup) != 0){
    error("Failed to read electric field data from file: %s\n", 
	  setup->field_name);
    return -1;
  }
  if (setup_wp(setup) != 0){
    error("Failed to read weighting potential from file %s\n",
	  setup->wp_name);
    return -1;
  }

  return 0;
}

static int efield_exists(cyl_pt pt, MJD_Siggen_Setup *setup){
  cyl_int_pt ipt;
  char ptstr[MAX_LINE];
  int  i, j, ir, iz;
  sprintf(ptstr, "(r,z) = (%.1f,%.1f)", pt.r, pt.z);
  if (outside_detector_cyl(pt, setup)){
    TELL_CHATTY("point %s is outside crystal\n", ptstr);
    return 0;
  }
  ipt.r = (pt.r - setup->rmin)/setup->rstep;  // CHECKED: no need for lrintf
  ipt.phi = 0;
  ipt.z = (pt.z - setup->zmin)/setup->zstep;  // CHECKED: no need for lrintf

  if (ipt.r < 0 || ipt.r + 1 >= setup->rlen ||
      ipt.z < 0 || ipt.z + 1 >= setup->zlen){
    TELL_CHATTY("point %s is outside wp table\n", ptstr);
    return 0;
  }
  for (i = 0; i < 2 ; i++){
    ir = ipt.r + i;
    for (j = 0; j < 2; j++){
      iz = ipt.z + j;
      if (setup->efld[ir][iz].r == 0.0 && setup->efld[ir][iz].z == 0.0) {
	TELL_CHATTY("point %s has no efield\n", ptstr);
	return 0;
      }
    }
  }
  TELL_CHATTY("point %s is in crystal\n", ptstr);
  return 1;
}

/* wpotential
   gives (interpolated) weighting potential at point pt, stored in wp
   returns 0 for success, 1 on failure
*/
int wpotential(point pt, float *wp, MJD_Siggen_Setup *setup){
  float w[2][2];
  int   i, j;
  cyl_int_pt ipt;
  cyl_pt cyl;

  // cyl = cart_to_cyl(pt);  // do not need to know phi, so save call to atan
  cyl.r = sqrt(pt.x*pt.x + pt.y*pt.y);
  cyl.z = pt.z;

  if (nearest_field_grid_index(cyl, &ipt, setup) < 0) return 1;
  grid_weights(cyl, ipt, w, setup);
  *wp = 0.0;
  for (i = 0; i < 2; i++){
    for (j = 0; j < 2; j++){
      *wp += w[i][j]*setup->wpot[ipt.r+i][ipt.z+j];
    }
  }

  return 0;
}

/* drift_velocity
   calculates drift velocity for charge q at point pt
   returns 0 on success, 1 on success but extrapolation was necessary,
   and -1 for failure
   anisotropic drift: crystal axes are assumed to be (x,y,z)
*/
int drift_velocity(point pt, float q, vector *velo, MJD_Siggen_Setup *setup){
  point  cart_en;
  cyl_pt e, en, cyl;
  cyl_int_pt ipt;
  int   i, sign;
  float abse, absv, f, a, b, c;
  float bp, cp, en4, en6;
  struct velocity_lookup *v_lookup1, *v_lookup2;

  /*  DCR: replaced this with faster code below, saves calls to atan and tan
  cyl = cart_to_cyl(pt);
  if (nearest_field_grid_index(cyl, &ipt, setup) < 0) return -1;
  e = efield(cyl, ipt, setup);
  abse = vector_norm_cyl(e, &en);
  en.phi = cyl.phi;
  cart_en = cyl_to_cart(en);
  */
  cyl.r = sqrt(pt.x*pt.x + pt.y*pt.y);
  cyl.z = pt.z;
  cyl.phi = 0;
  if (nearest_field_grid_index(cyl, &ipt, setup) < 0) return -1;
  e = efield(cyl, ipt, setup);
  abse = vector_norm_cyl(e, &en);
  if (cyl.r > 0.001) {
    cart_en.x = en.r * pt.x/cyl.r;
    cart_en.y = en.r * pt.y/cyl.r;
  } else {
    cart_en.x = cart_en.y = 0;
  }
  cart_en.z = en.z;

  /* find location in table to interpolate from */
  for (i = 0; i < setup->v_lookup_len - 2 && abse > setup->v_lookup[i+1].e; i++);
  v_lookup1 = setup->v_lookup + i;
  v_lookup2 = setup->v_lookup + i+1;
  f = (abse - v_lookup1->e)/(v_lookup2->e - v_lookup1->e);
  if (q > 0){
    a = (v_lookup2->ha - v_lookup1->ha)*f+v_lookup1->ha;
    b = (v_lookup2->hb- v_lookup1->hb)*f+v_lookup1->hb;
    c = (v_lookup2->hc - v_lookup1->hc)*f+v_lookup1->hc;
    bp = (v_lookup2->hbp- v_lookup1->hbp)*f+v_lookup1->hbp;
    cp = (v_lookup2->hcp - v_lookup1->hcp)*f+v_lookup1->hcp;
    setup->dv_dE = (v_lookup2->h100 - v_lookup1->h100)/(v_lookup2->e - v_lookup1->e);
  }else{
    a = (v_lookup2->ea - v_lookup1->ea)*f+v_lookup1->ea;
    b = (v_lookup2->eb- v_lookup1->eb)*f+v_lookup1->eb;
    c = (v_lookup2->ec - v_lookup1->ec)*f+v_lookup1->ec;
    bp = (v_lookup2->ebp- v_lookup1->ebp)*f+v_lookup1->ebp;
    cp = (v_lookup2->ecp - v_lookup1->ecp)*f+v_lookup1->ecp;
    setup->dv_dE = (v_lookup2->e100 - v_lookup1->e100)/(v_lookup2->e - v_lookup1->e);
  }
  /* velocity can vary from the direction of the el. field
     due to effect of crystal axes */
#define POW4(x) ((x)*(x)*(x)*(x))
#define POW6(x) ((x)*(x)*(x)*(x)*(x)*(x))
  en4 = POW4(cart_en.x) + POW4(cart_en.y) + POW4(cart_en.z);
  en6 = POW6(cart_en.x) + POW6(cart_en.y) + POW6(cart_en.z);
  absv = a + b*en4 + c*en6;
  sign = (q < 0 ? -1 : 1);
  setup->v_over_E = absv / abse;
  velo->x = sign*cart_en.x*(absv+bp*4*(cart_en.x*cart_en.x - en4)
			    + cp*6*(POW4(cart_en.x) - en6));
  velo->y = sign*cart_en.y*(absv+bp*4*(cart_en.y*cart_en.y - en4)
			    + cp*6*(POW4(cart_en.y) - en6));
  velo->z = sign*cart_en.z*(absv+bp*4*(cart_en.z*cart_en.z - en4)
			    + cp*6*(POW4(cart_en.z) - en6));
#undef POW4
#undef POW6
  return 0;
}

/* Find (interpolated or extrapolated) electric field for this point */
static cyl_pt efield(cyl_pt pt, cyl_int_pt ipt, MJD_Siggen_Setup *setup){
  cyl_pt e = {0,0,0}, ef;
  float  w[2][2];
  int    i, j;

  grid_weights(pt, ipt, w, setup);
  for (i = 0; i < 2; i++){
    for (j = 0; j < 2; j++){
      ef = setup->efld[ipt.r + i][ipt.z + j];
      e.r += ef.r*w[i][j];
      e.z += ef.z*w[i][j];
    }
  }
  e.phi = pt.phi;
  return e;
}


/* Find weights for 8 voxel corner points around pt for e/wp field*/
/* DCR: modified to work for both interpolation and extrapolation */
static int grid_weights(cyl_pt pt, cyl_int_pt ipt, float out[2][2],
			MJD_Siggen_Setup *setup){
  float r, z;

  r = (pt.r - setup->rmin)/setup->rstep - ipt.r;
  z = (pt.z - setup->zmin)/setup->zstep - ipt.z;

  out[0][0] = (1.0 - r) * (1.0 - z);
  out[0][1] = (1.0 - r) *        z;
  out[1][0] =        r  * (1.0 - z);
  out[1][1] =        r  *        z;
  return 0;
}


/*find existing integer field grid index closest to pt*/
/* added DCR */
static int nearest_field_grid_index(cyl_pt pt, cyl_int_pt *ipt,
				    MJD_Siggen_Setup *setup){
  /* returns <0 if outside crystal or too far from a valid grid point
              0 if interpolation is okay
              1 if we can find a point but extrapolation is needed
  */
  static cyl_pt  last_pt;
  static cyl_int_pt last_ipt;
  static int     last_ret = -99;
  cyl_pt new_pt;
  int    dr, dz;
  float  d[3] = {0.0, -1.0, 1.0};

  if (last_ret != -99 &&
      pt.r == last_pt.r && pt.z == last_pt.z) {
    *ipt = last_ipt;
    return last_ret;
  }
  last_pt = pt;
  last_ret = -2;

  if (outside_detector_cyl(pt, setup)) {
    last_ret = -1;
  } else{
    new_pt.phi = 0.0;
    for (dz=0; dz<3; dz++) {
      new_pt.z = pt.z + d[dz]*setup->zstep;
      for (dr=0; dr<3; dr++) {
	new_pt.r = pt.r + d[dr]*setup->rstep;
	if (efield_exists(new_pt, setup)) {
	  last_ipt.r = (new_pt.r - setup->rmin)/setup->rstep;  // CHECKED: do NOT use lrintf
	  last_ipt.phi = 0;
	  last_ipt.z = (new_pt.z - setup->zmin)/setup->zstep;  // CHECKED: do NOT use lrintf
	  *ipt = last_ipt;
	  if (dr == 0 && dz == 0) {
	    last_ret = 0;
	  } else {
	    last_ret = 1;
	  }
	  return last_ret;
	}
      }
    }
  }

  return last_ret;
}

/* setup_velo
   set up drift velocity calculations (read in table)
*/
static int setup_velo(MJD_Siggen_Setup *setup){
  static int vlook_sz = 0;
  static struct velocity_lookup *v_lookup;

  char  line[MAX_LINE], *c;
  FILE  *fp;
  int   i, v_lookup_len;
  struct velocity_lookup *tmp, v, v0;
  float sumb_e, sumc_e, sumb_h, sumc_h;

  double be=1.3e7, bh=1.2e7, thetae=200.0, thetah=200.0;  // parameters for temperature correction
  double pwre=-1.680, pwrh=-2.398, mue=5.66e7, muh=1.63e9; //     adopted for Ge   DCR Feb 2015
  double mu_0_1, mu_0_2, v_s_1, v_s_2, E_c_1, E_c_2, e, f;

  if (vlook_sz == 0) {
    vlook_sz = 10;
    if ((v_lookup = (struct velocity_lookup *)
	 malloc(vlook_sz*sizeof(*v_lookup))) == NULL) {
      error("malloc failed in setup_velo\n");
      return -1;
    }
  }
  if ((fp = fopen(setup->drift_name, "r")) == NULL){
    error("failed to open velocity lookup table file: '%s'\n", setup->drift_name);
    return -1;
  }
  line[0] = '#';
  c = line;
  while ((line[0] == '#' || line[0] == '\0') && c != NULL) c = fgets(line, MAX_LINE, fp);
  if (c == NULL) {
    error("Failed to read velocity lookup table from file: %s\n", setup->drift_name);
    fclose(fp);
    return -1;
  }
  TELL_CHATTY("Drift velocity table:\n"
	      "  e          e100    e110    e111    h100    h110    h111\n");   
  for (v_lookup_len = 0; ;v_lookup_len++){
    if (v_lookup_len == vlook_sz - 1){
      vlook_sz += 10;
      if ((tmp = (struct velocity_lookup *)
	   realloc(v_lookup, vlook_sz*sizeof(*v_lookup))) == NULL){
	error("realloc failed in setup_velo\n");
	fclose(fp);
	return -1;
      }
      v_lookup = tmp;
    }
    if (sscanf(line, "%f %f %f %f %f %f %f", 
	       &v_lookup[v_lookup_len].e,
	       &v_lookup[v_lookup_len].e100,
	       &v_lookup[v_lookup_len].e110,
	       &v_lookup[v_lookup_len].e111,
	       &v_lookup[v_lookup_len].h100,
	       &v_lookup[v_lookup_len].h110,
	       &v_lookup[v_lookup_len].h111) != 7){
      break; //assume EOF
    }	   
    //v_lookup[v_lookup_len].e *= 100; /*V/m*/
    tmp = &v_lookup[v_lookup_len];
    TELL_CHATTY("%10.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f\n",
		tmp->e, tmp->e100, tmp->e110, tmp->e111, tmp->h100, tmp->h110,tmp->h111);
    line[0] = '#';
    while ((line[0] == '#' || line[0] == '\0' ||
	    line[0] == '\n' || line[0] == '\r') && c != NULL) c = fgets(line, MAX_LINE, fp);
    if (c == NULL) break;
    if (line[0] == 'e' || line[0] == 'h') break; /* no more velocities data;
						    now reading temp correction data */
  }

  /* check for and decode temperature correction parameters */
  while (line[0] == 'e' || line[0] == 'h') {
    if (line[0] == 'e' &&
	sscanf(line+2, "%lf %lf %lf %lf", 
	       &mue, &pwre, &be, &thetae) != 4) break;//asume EOF
    if (line[0] == 'h' &&
	sscanf(line+2, "%lf %lf %lf %lf", 
	       &muh, &pwrh, &bh, &thetah) != 4) break;//asume EOF
    if (line[0] == 'e')
      TELL_CHATTY("electrons: mu_0 = %.2e x T^%.4f  B = %.2e  Theta = %.0f\n",
		  mue, pwre, be, thetae);
    if (line[0] == 'h')
      TELL_CHATTY("    holes: mu_0 = %.2e x T^%.4f  B = %.2e  Theta = %.0f\n",
		  muh, pwrh, bh, thetah);

    line[0] = '#';
    while ((line[0] == '#' || line[0] == '\0') && c != NULL) c = fgets(line, MAX_LINE, fp);
    if (c == NULL) break;
  }

  if (v_lookup_len == 0){
    error("Failed to read velocity lookup table from file: %s\n", setup->drift_name);
    return -1;
  }  
  v_lookup_len++;
  if (vlook_sz != v_lookup_len){
    if ((tmp = (struct velocity_lookup *) 
	 realloc(v_lookup, v_lookup_len*sizeof(*v_lookup))) == NULL){
      error("realloc failed in setup_velo. This should not happen\n");
      fclose(fp);
      return -1;
    }
    v_lookup = tmp;
    vlook_sz = v_lookup_len;
  }
  TELL_NORMAL("Drift velocity table has %d rows of data\n", v_lookup_len);
  fclose(fp);

  /*
    apply temperature dependence to mobilities;
    see drift_velocities.doc and tempdep.c
    The drift velocity reduces at higher temperature due to the increasing of
    scattering with the lattice vibration. We used a model by M. Ali Omar and
    L. Reggiani (Solid-State Electronics Vol. 30, No. 12 (1987) 1351) to
    calculate the temperature dependence.
  */
  /* electrons */
  TELL_NORMAL("Adjusting mobilities for temperature, from %.1f to %.1f\n", REF_TEMP, setup->xtal_temp);
  TELL_CHATTY("Index  field  vel_factor\n");
  mu_0_1 = mue * pow(REF_TEMP, pwre);
  v_s_1 = be * sqrt(tanh(0.5 * thetae / REF_TEMP));
  E_c_1 = v_s_1 / mu_0_1;
  mu_0_2 = mue * pow(setup->xtal_temp, pwre);
  v_s_2 = be * sqrt(tanh(0.5 * thetae / setup->xtal_temp));
  E_c_2 = v_s_2 / mu_0_2;
  for (i = 0; i < vlook_sz; i++){
    e = v_lookup[i].e;
    if (e < 1) continue;
    f = (v_s_2 * (e/E_c_2) / sqrt(1.0 + (e/E_c_2) * (e/E_c_2))) /
        (v_s_1 * (e/E_c_1) / sqrt(1.0 + (e/E_c_1) * (e/E_c_1)));
    v_lookup[i].e100 *= f;
    v_lookup[i].e110 *= f;
    v_lookup[i].e111 *= f;
    TELL_CHATTY("%2d %5.0f %f\n", i, e, f);
  }

  /* holes */
  mu_0_1 = muh * pow(REF_TEMP, pwrh);
  v_s_1 = bh * sqrt(tanh(0.5 * thetah / REF_TEMP));
  E_c_1 = v_s_1 / mu_0_1;
  mu_0_2 = muh * pow(setup->xtal_temp, pwrh);
  v_s_2 = bh * sqrt(tanh(0.5 * thetah / setup->xtal_temp));
  E_c_2 = v_s_2 / mu_0_2;
  for (i = 0; i < vlook_sz; i++){
    e = v_lookup[i].e;
    if (e < 1) continue;
    f = (v_s_2 * (e/E_c_2) / sqrt(1.0 + (e/E_c_2) * (e/E_c_2))) /
        (v_s_1 * (e/E_c_1) / sqrt(1.0 + (e/E_c_1) * (e/E_c_1)));
    v_lookup[i].h100 *= f;
    v_lookup[i].h110 *= f;
    v_lookup[i].h111 *= f;
    TELL_CHATTY("%2d %5.0f %f\n", i, e, f);
  }
  /* end of temperature correction */

  for (i = 0; i < vlook_sz; i++){
    v = v_lookup[i];
    v_lookup[i].ea =  0.5 * v.e100 -  4 * v.e110 +  4.5 * v.e111;
    v_lookup[i].eb = -2.5 * v.e100 + 16 * v.e110 - 13.5 * v.e111;
    v_lookup[i].ec =  3.0 * v.e100 - 12 * v.e110 +  9.0 * v.e111;
    v_lookup[i].ha =  0.5 * v.h100 -  4 * v.h110 +  4.5 * v.h111;
    v_lookup[i].hb = -2.5 * v.h100 + 16 * v.h110 - 13.5 * v.h111;
    v_lookup[i].hc =  3.0 * v.h100 - 12 * v.h110 +  9.0 * v.h111;
  }
  v_lookup[0].ebp = v_lookup[0].ecp = v_lookup[0].hbp = v_lookup[0].hcp = 0.0;
  sumb_e = sumc_e = sumb_h = sumc_h = 0.0;
  for (i = 1; i < vlook_sz; i++){
    v0 = v_lookup[i-1];
    v = v_lookup[i];
    sumb_e += (v.e - v0.e)*(v0.eb+v.eb)/2;
    sumc_e += (v.e - v0.e)*(v0.ec+v.ec)/2;
    sumb_h += (v.e - v0.e)*(v0.hb+v.hb)/2;
    sumc_h += (v.e - v0.e)*(v0.hc+v.hc)/2;
    v_lookup[i].ebp = sumb_e/v.e;
    v_lookup[i].ecp = sumc_e/v.e;
    v_lookup[i].hbp = sumb_h/v.e;
    v_lookup[i].hcp = sumc_h/v.e;
  }

  setup->v_lookup = v_lookup;
  setup->v_lookup_len = v_lookup_len;

  return 0;
}


/* This may or may not break if we switch to a non-integer grid*/
/*setup_efield
  read electric field data from file, apply sanity checks
  returns 0 for success
*/
static int setup_efield(MJD_Siggen_Setup *setup){
  FILE   *fp;
  char   line[MAX_LINE], *cp;
  int    i, j, lineno = 0;
  float  v, eabs, er, ez;
  cyl_pt cyl, **efld;

  if ((fp = fopen(setup->field_name, "r")) == NULL){
    error("failed to open electric field table: %s\n", setup->field_name);
    return 1;
  }
  
  setup->rlen = lrintf((setup->rmax - setup->rmin)/setup->rstep) + 1;
  setup->zlen = lrintf((setup->zmax - setup->zmin)/setup->zstep) + 1;
  TELL_CHATTY("rlen, zlen: %d, %d\n", setup->rlen, setup->zlen);

  // here I assume that r, zlen never change from their initial values, which is reasonable
  if ((efld = (cyl_pt **) malloc(setup->rlen*sizeof(*efld))) == NULL) {
    error("Malloc failed in setup_efield\n");
    fclose(fp);
    return 1;
  }
  for (i = 0; i < setup->rlen; i++){
    if ((efld[i] = (cyl_pt *) malloc(setup->zlen*sizeof(*efld[i]))) == NULL){
      error("Malloc failed in setup_efield\n");
      //NB: potential memory leak here.
      fclose(fp);
      return 1;
    }
    memset(efld[i], 0, setup->zlen*sizeof(*efld[i]));
  }
  TELL_NORMAL("Reading electric field data from file: %s\n", setup->field_name);

  /*now read the table*/
  while(fgets(line, MAX_LINE, fp) != NULL){
    lineno++;
    for (cp = line; isspace(*cp) && *cp != '\0'; cp++);
    if (*cp == '#' || !strlen(cp)) continue;
    if (sscanf(line, "%f %f %f %f %f %f", 
	       &cyl.r, &cyl.z, &v, &eabs, &er, &ez) != 6){
      error("failed to read electric field data from line no %d\n"
	    "of file %s\n", lineno, setup->field_name);
      fclose(fp);
      return 1;
    }
    i = lrintf((cyl.r - setup->rmin)/setup->rstep);
    j = lrintf((cyl.z - setup->zmin)/setup->zstep);
    if (i < 0 || i >= setup->rlen || j < 0 || j >= setup->zlen) {
      // error("Error in efield line %d, i = %d, j = %d\n", line, i, j);
      continue;
    }
    cyl.phi = 0;
    if (outside_detector_cyl(cyl, setup)) continue;
    efld[i][j].r = er;
    efld[i][j].z = ez;
    efld[i][j].phi = 0;
  }      

  TELL_NORMAL("Done reading %d lines of electric field data\n", lineno);
  fclose(fp);

  setup->efld = efld;
  for (i = 0; i < setup->rlen; i++) setup->efld[i] = efld[i];

  return 0;
}

/*setup_wp
  read weighting potential values from files. returns 0 on success*/
static int setup_wp(MJD_Siggen_Setup *setup){
  FILE   *fp;
  char   line[MAX_LINE], *cp;
  int    i, j, lineno;
  cyl_pt cyl;
  float  wp, **wpot;

  setup->rlen = lrintf((setup->rmax - setup->rmin)/setup->rstep) + 1;
  setup->zlen = lrintf((setup->zmax - setup->zmin)/setup->zstep) + 1;
  TELL_CHATTY("rlen, zlen: %d, %d\n", setup->rlen, setup->zlen);

  //assuming rlen, zlen never change as for setup_efld
  if ((wpot = (float **) malloc(setup->rlen*sizeof(*wpot))) == NULL){
    error("Malloc failed in setup_wp\n");
    return 1;
  }
  for (i = 0; i < setup->rlen; i++){
    if ((wpot[i] = (float *) malloc(setup->zlen*sizeof(*wpot[i]))) == NULL){  
      error("Malloc failed in setup_wp\n");
      //NB: memory leak here.
      return 1;
    }
    memset(wpot[i], 0, setup->zlen*sizeof(*wpot[i]));
  }
  if ((fp = fopen(setup->wp_name, "r")) == NULL){
    error("failed to open file: %s\n", setup->wp_name);
    return -1;
  }
  lineno = 0;
  TELL_NORMAL("Reading weighting potential from file: %s\n", setup->wp_name);
  while (fgets(line, MAX_LINE, fp) != NULL){
    lineno++;
    for (cp = line; isspace(*cp) && *cp != '\0'; cp++);
    if (*cp == '#' || !strlen(cp)) continue;
    if (sscanf(line, "%f %f %f\n",&cyl.r, &cyl.z, &wp) != 3){ 
      error("failed to read weighting potential from line %d\n"
	    "line: %s", lineno, line);
      fclose(fp);
      return 1;
    }
    i = lrintf((cyl.r - setup->rmin)/setup->rstep);
    j = lrintf((cyl.z - setup->zmin)/setup->zstep);
    if (i < 0 || i >= setup->rlen || j < 0 || j >= setup->zlen) continue;
    if (outside_detector_cyl(cyl, setup)) continue;
    wpot[i][j] = wp;
  }
  TELL_NORMAL("Done reading %d lines of WP data\n", lineno);
  fclose(fp);

  setup->wpot = wpot;
  for (i = 0; i < setup->rlen; i++) setup->wpot[i] = wpot[i];

  return 0;
}


/* free malloc()'ed memory and do other cleanup*/
int fields_finalize(MJD_Siggen_Setup *setup){
  int i;

  for (i = 0; i < lrintf((setup->rmax - setup->rmin)/setup->rstep) + 1; i++){
    free(setup->efld[i]);
    free(setup->wpot[i]);
  }
  free(setup->efld);
  free(setup->wpot);
  free(setup->v_lookup);
  setup->efld = NULL;
  setup->wpot = NULL;
  setup->v_lookup = NULL;

  return 1;
}

void set_temp(float temp, MJD_Siggen_Setup *setup){
  if (temp < MIN_TEMP || temp > MAX_TEMP){
    error("temperature out of range: %f\n", temp);
  }else{
    setup->xtal_temp = temp;
    error("temperature set to %f\n", temp);
    /* re-read velocities and correct them to the new temperature value */
    setup_velo(setup);
  }
}
