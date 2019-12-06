/* fields.c -- based on m3d2s.f by I-Yang Lee
 * Karin Lagergren
 *
 * This module handles the electric field and weighting potential and 
 * calculates drift velocities
 *
 * November 2007 -- setup_efield, setup_wp rewritten
 */

//TODO: Add setup_done flag & check it before trying to access data ?
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "point.h"
#include "fields.h"
#include "detector_geometry.h"
#include "calc_signal.h"
#include "signal_calc_util.h"

static int    nearest_field_grid_index(MJD_Siggen_Setup *setup, point pt, int_pt *ipt);
static int    grid_weights(MJD_Siggen_Setup *setup, point pt, int_pt ipt, float out[2][2][2]);
static vector efield(MJD_Siggen_Setup *setup, point pt, int_pt ipt);
static int    setup_efield(MJD_Siggen_Setup *setup);
static int    setup_wp(MJD_Siggen_Setup *setup);
static int    setup_velo(MJD_Siggen_Setup *setup);
static int_pt field_grid_index(MJD_Siggen_Setup *setup, point pt);


/* field_setup
   given a field directory file, read electic field and weighting
   potential tables from files listed in directory
   returns 0 for success
*/
int field_setup(MJD_Siggen_Setup *setup) {

  setup->xmin = setup->ymin = -setup->xtal_radius;
  setup->xmax = setup->ymax = setup->xtal_radius;
  setup->zmin = 0; setup->zmax = setup->xtal_length;

  setup->numx = (int)rint((setup->xmax - setup->xmin)/setup->xtal_grid) + 1;
  setup->numy = (int)rint((setup->ymax - setup->ymin)/setup->xtal_grid) + 1;
  setup->numz = (int)rint((setup->zmax - setup->zmin)/setup->xtal_grid) + 1;
  tell(NORMAL, "numx, numy, numz: %d %d %d\n", setup->numx, setup->numy, setup->numz);

  tell(NORMAL, "Detector temperature is set to %.1f K\n", setup->xtal_temp);

  if (setup_velo(setup) != 0) {
    error("failed to read drift velocity data from file: %s\n", setup->drift_name);
    return -1;
  }

  if (setup_efield(setup) != 0) {
    error("Failed to read electric field data from file: %s\n",
	  setup->field_name);
    return -1;
  }
  tell(NORMAL, "===> OK\n");    
  if (setup_wp(setup) != 0) {
    error("failed to set up weighting potentials\n");
    return -1;
  }
  return 0;
}


static int efield_exists(MJD_Siggen_Setup *setup, point pt) {
  char ptstr[MAX_LINE];
  int i, j, k;
  int ix, iy, iz;
  int_pt ipt;

  pt_to_str(ptstr, MAX_LINE, pt);
  if (!in_crystal(setup, pt)) {
    tell(CHATTY, "point %s is outside crystal\n", ptstr);
    return 0;
  } else if (0) {
    tell(CHATTY, "point %s is in crystal\n", ptstr);
  }
  ipt = field_grid_index(setup, pt);

  if (ipt.x < 0 || ipt.x + 1 >= setup->numx ||
      ipt.y < 0 || ipt.y + 1 >= setup->numy ||
      ipt.z < 0 || ipt.z + 1 >= setup->numz) {
    tell(CHATTY, "point %s is outside wp table\n", ptstr);
    return 0;
  }
  for (i = 0; i < 2 ; i++) {
    ix = ipt.x + i;
    for (j = 0; j < 2; j++) {
      iy = ipt.y + j;
      for (k = 0; k < 2; k++) {
	iz = ipt.z + k;
	if (setup->efld[ix][iy][iz].x == 0.0 && setup->efld[ix][iy][iz].y == 0.0 
	    && setup->efld[ix][iy][iz].z == 0.0 ) return 0;
      }
    }
  }

  return 1;
}


/* wpotential
   gives (interpolated) weighting potential 
   at point pt. returns 0 for success, 1 on failure
*/
int wpotential(MJD_Siggen_Setup *setup, point pt, float *wp) {
  float w[2][2][2];
  int i, j, k;
  int_pt ipt;
  int res;

  res = nearest_field_grid_index(setup, pt, &ipt);
  if (res < 0) return 1;
  grid_weights(setup, pt, ipt, w);
  *wp = 0.0;
  for (i = 0; i < 2; i++) {
    for (j = 0; j < 2; j++) {
      for (k = 0; k < 2; k++) {
        *wp += w[i][j][k] * setup->wpot[ipt.x+i][ipt.y+j][ipt.z+k];
      }
    }
  }
  return 0;
}

/* drift_velocity
   calculates drift velocity for charge q at point pt
   returns 0 on success, > 0 on success but extrapolation was necessary,
   and -1 for failure
   anisotropic drift: crystal axes are assumed to be (x,y,z)
*/
int drift_velocity(MJD_Siggen_Setup *setup, point pt, float q, vector *velo) {
  vector e, en;
  int_pt ipt;
  int    i;
  float  abse;
  float  f;
  float  a, b, c;
  float  absv;
  int    sign;
  float  bp, cp;
  float  en4, en6;
  vector v;
  int    nfgi_ans;

  if ((nfgi_ans = nearest_field_grid_index(setup, pt, &ipt)) < 0) return -1;
  /*   nfgi_ans < 0 if outside crystal or too far from a valid grid point
                = 0 if interpolation is okay
                > 1 if we can find a point but extrapolation is needed */
  tell(CHATTY, "field grid exists for point %d %d %d -- (%.2f %.2f %.2f)\n",
       ipt.x, ipt.y, ipt.z, pt.x, pt.y, pt.z);
  if (nfgi_ans)
    tell(CHATTY, "  ... but need to extrapolate.\n");
  e = efield(setup, pt, ipt);
  abse = vector_norm(e, &en);
 
  /* find location in table to interpolate / extrapolate from*/
  for (i = 0; i < setup->v_lookup_len - 2 && abse > setup->v_lookup[i+1].e; i++);
  /* interpolate / extrapolate */
  f = (abse - setup->v_lookup[i].e)/(setup->v_lookup[i+1].e - setup->v_lookup[i].e); 
  if (q > 0) {  /* hole */
    a = (setup->v_lookup[i+1].ha - setup->v_lookup[i].ha)*f+setup->v_lookup[i].ha;
    b = (setup->v_lookup[i+1].hb- setup->v_lookup[i].hb)*f+setup->v_lookup[i].hb;
    c = (setup->v_lookup[i+1].hc - setup->v_lookup[i].hc)*f+setup->v_lookup[i].hc;
    bp = (setup->v_lookup[i+1].hbp- setup->v_lookup[i].hbp)*f+setup->v_lookup[i].hbp;
    cp = (setup->v_lookup[i+1].hcp - setup->v_lookup[i].hcp)*f+setup->v_lookup[i].hcp;
  } else {   /* electron */
    a = (setup->v_lookup[i+1].ea - setup->v_lookup[i].ea)*f+setup->v_lookup[i].ea;
    b = (setup->v_lookup[i+1].eb- setup->v_lookup[i].eb)*f+setup->v_lookup[i].eb;
    c = (setup->v_lookup[i+1].ec - setup->v_lookup[i].ec)*f+setup->v_lookup[i].ec;
    bp = (setup->v_lookup[i+1].ebp- setup->v_lookup[i].ebp)*f+setup->v_lookup[i].ebp;
    cp = (setup->v_lookup[i+1].ecp - setup->v_lookup[i].ecp)*f+setup->v_lookup[i].ecp;
  }

#define POW4(x) ((x)*(x)*(x)*(x))
#define POW6(x) ((x)*(x)*(x)*(x)*(x)*(x))
  en4 = POW4(en.x) + POW4(en.y) + POW4(en.z);
  en6 = POW6(en.x) + POW6(en.y) + POW6(en.z);
  absv = a + b*en4 + c*en6;
  sign = (q < 0 ? -1 : 1);
  velo->x = sign*en.x*(absv+bp*4*(en.x*en.x - en4) + cp*6*(POW4(en.x) - en6));
  velo->y = sign*en.y*(absv+bp*4*(en.y*en.y - en4) + cp*6*(POW4(en.y) - en6));
  velo->z = sign*en.z*(absv+bp*4*(en.z*en.z - en4) + cp*6*(POW4(en.z) - en6));
  v = *velo;
  velo->x = v.x;
  velo->y = v.y;
  velo->z = v.z;
  
#undef POW4
#undef POW6
  
  return nfgi_ans;
}

/* Find (interpolated or extrapolated) electric field for this point */
static vector efield(MJD_Siggen_Setup *setup, point pt, int_pt ipt) {
  vector e, zero = {0,0,0}, ef;
  float w[2][2][2];
  int i, j, k;

  grid_weights(setup, pt, ipt, w);
  e = zero;
  for (i = 0; i < 2; i++) {
    for (j = 0; j < 2; j++) {
      for (k = 0; k < 2; k++) {
	ef = setup->efld[ipt.x + i][ipt.y + j][ipt.z +k];
	e.x += ef.x*w[i][j][k];
	e.y += ef.y*w[i][j][k];
	e.z += ef.z*w[i][j][k];
      }
    }
  }
  
  return e;
}


/* Find weights for 8 voxel corner points around pt for e/wp field*/
/* DCR: modified to work for both interpolation and extrapolation */
static int grid_weights(MJD_Siggen_Setup *setup, point pt, int_pt ipt, float out[2][2][2]) {
  float x, y, z;

  x = (pt.x - setup->xmin)/setup->xtal_grid - ipt.x;
  y = (pt.y - setup->ymin)/setup->xtal_grid - ipt.y;
  z = (pt.z - setup->zmin)/setup->xtal_grid - ipt.z;

  out[0][0][0] = (1.0 - x) * (1.0 - y) * (1.0 - z);
  out[0][0][1] = (1.0 - x) * (1.0 - y) *        z;
  out[0][1][0] = (1.0 - x) *        y  * (1.0 - z);
  out[0][1][1] = (1.0 - x) *        y  *        z;
  out[1][0][0] =        x  * (1.0 - y) * (1.0 - z);
  out[1][0][1] =        x  * (1.0 - y) *        z;
  out[1][1][0] =        x  *        y  * (1.0 - z);
  out[1][1][1] =        x  *        y  *        z;

  return 0;
}


/*find existing integer field grid index closest to pt*/
/* added DCR */
static int nearest_field_grid_index(MJD_Siggen_Setup *setup, point pt, int_pt *ipt) {
  /* returns < 0 if outside crystal or too far from a valid grid point
               0 if interpolation is okay
             > 0 if we can find a point but extrapolation is needed
  */
  static point  last_pt;
  static int_pt last_ipt;
  static int    last_ret = -99;
  point new;
  int x, y, z;
  float d[3] = {0.0, -1.0, 1.0};

  if (last_ret != -99 &&
      pt.x == last_pt.x && pt.y == last_pt.y && pt.z == last_pt.z) {
    *ipt = last_ipt;
    return last_ret;
  }
  last_pt = pt;
  last_ret = -2;

  if (!in_crystal(setup, pt)) {
    last_ret = -1;
  } else if (fabs(pt.x) > fabs(pt.y)) {
    for (z=0; z<3; z++) {
      new.z = pt.z + d[z] * setup->xtal_grid;
      for (y=0; y<3; y++) {
	new.y = pt.y + d[y] * setup->xtal_grid;
	for (x=0; x<3; x++) {
	  new.x = pt.x + d[x] * setup->xtal_grid;
	  if (efield_exists(setup, new)) {
	    *ipt = last_ipt = field_grid_index(setup, new);
	    if (x == 0 && y == 0 && z == 0) {
	      last_ret = 0;
	    } else {
	      last_ret = 1;
	    }
	    return last_ret;
	  }
	}
      }
    }

  } else {
    for (z=0; z<3; z++) {
      new.z = pt.z + d[z] * setup->xtal_grid;
      for (x=0; x<3; x++) {
	new.x = pt.x + d[x] * setup->xtal_grid;
	for (y=0; y<3; y++) {
	  new.y = pt.y + d[y] * setup->xtal_grid;
	  if (efield_exists(setup, new)) {
	    *ipt = last_ipt = field_grid_index(setup, new);
	    if (x == 0 && y == 0 && z == 0) {
	      last_ret = 0;
	    } else {
	      last_ret = 1;
	    }
	    return last_ret;
	  }
	}
      }
    }
  }

  return last_ret;
}

/*find integer field grid index corresponding to pt*/
static int_pt field_grid_index(MJD_Siggen_Setup *setup, point pt) {
  int_pt ipt;

  ipt.x = (pt.x - setup->xmin)/setup->xtal_grid;  // NOTE: tried adding lrintf()  Oct2019
  ipt.y = (pt.y - setup->ymin)/setup->xtal_grid;  // CHECKED - lrintf causes signal calc failures
  ipt.z = (pt.z - setup->zmin)/setup->xtal_grid;  //   and bad field. Do NOT add lrintf here.

  return ipt;
}

/* setup_velo
   set up drift velocity calculations (read in table)
*/
static int setup_velo(MJD_Siggen_Setup *setup) {
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
  if ((fp = fopen(setup->drift_name, "r")) == NULL) {
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
  tell(CHATTY, "Drift velocity table:\n"
	       "  e          e100    e110    e111    h100    h110    h111\n");   
  for (v_lookup_len = 0; ;v_lookup_len++) {
    if (v_lookup_len == vlook_sz - 1) {
      vlook_sz += 10;
      if ((tmp = (struct velocity_lookup *)
	   realloc(v_lookup, vlook_sz*sizeof(*v_lookup))) == NULL) {
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
	       &v_lookup[v_lookup_len].h111) != 7) {
      break; //assume EOF
    }	   
    //v_lookup[v_lookup_len].e *= 100; /*V/m*/
    tmp = &v_lookup[v_lookup_len];
    tell(CHATTY, "%10.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f\n",
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
      tell(CHATTY, "electrons: mu_0 = %.2e x T^%.4f  B = %.2e  Theta = %.0f\n",
           mue, pwre, be, thetae);
    if (line[0] == 'h')
      tell(CHATTY, "    holes: mu_0 = %.2e x T^%.4f  B = %.2e  Theta = %.0f\n",
           muh, pwrh, bh, thetah);

    line[0] = '#';
    while ((line[0] == '#' || line[0] == '\0') && c != NULL) c = fgets(line, MAX_LINE, fp);
    if (c == NULL) break;
  }

  if (v_lookup_len == 0) {
    error("Failed to read velocity lookup table from file: %s\n", setup->drift_name);
    return -1;
  }  
  v_lookup_len++;
  if (vlook_sz != v_lookup_len) {
    if ((tmp = (struct velocity_lookup *) 
	 realloc(v_lookup, v_lookup_len*sizeof(*v_lookup))) == NULL) {
      error("realloc failed in setup_velo. This should not happen\n");
      fclose(fp);
      return -1;
    }
    v_lookup = tmp;
    vlook_sz = v_lookup_len;
  }
  tell(NORMAL, "Drift velocity table has %d rows of data\n", v_lookup_len);
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
  tell(NORMAL, "Adjusting mobilities for temperature, from %.1f to %.1f\n", REF_TEMP, setup->xtal_temp);
  tell(CHATTY, "Index  field  vel_factor\n");
  mu_0_1 = mue * pow(REF_TEMP, pwre);
  v_s_1 = be * sqrt(tanh(0.5 * thetae / REF_TEMP));
  E_c_1 = v_s_1 / mu_0_1;
  mu_0_2 = mue * pow(setup->xtal_temp, pwre);
  v_s_2 = be * sqrt(tanh(0.5 * thetae / setup->xtal_temp));
  E_c_2 = v_s_2 / mu_0_2;
  for (i = 0; i < vlook_sz; i++) {
    e = v_lookup[i].e;
    if (e < 1) continue;
    f = (v_s_2 * (e/E_c_2) / sqrt(1.0 + (e/E_c_2) * (e/E_c_2))) /
        (v_s_1 * (e/E_c_1) / sqrt(1.0 + (e/E_c_1) * (e/E_c_1)));
    v_lookup[i].e100 *= f;
    v_lookup[i].e110 *= f;
    v_lookup[i].e111 *= f;
    tell(CHATTY, "%2d %5.0f %f\n", i, e, f);
  }

  /* holes */
  mu_0_1 = muh * pow(REF_TEMP, pwrh);
  v_s_1 = bh * sqrt(tanh(0.5 * thetah / REF_TEMP));
  E_c_1 = v_s_1 / mu_0_1;
  mu_0_2 = muh * pow(setup->xtal_temp, pwrh);
  v_s_2 = bh * sqrt(tanh(0.5 * thetah / setup->xtal_temp));
  E_c_2 = v_s_2 / mu_0_2;
  for (i = 0; i < vlook_sz; i++) {
    e = v_lookup[i].e;
    if (e < 1) continue;
    f = (v_s_2 * (e/E_c_2) / sqrt(1.0 + (e/E_c_2) * (e/E_c_2))) /
        (v_s_1 * (e/E_c_1) / sqrt(1.0 + (e/E_c_1) * (e/E_c_1)));
    v_lookup[i].h100 *= f;
    v_lookup[i].h110 *= f;
    v_lookup[i].h111 *= f;
    tell(CHATTY, "%2d %5.0f %f\n", i, e, f);
  }
  /* end of temperature correction */

  for (i = 0; i < vlook_sz; i++) {
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
  for (i = 1; i < vlook_sz; i++) {
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
static int setup_efield(MJD_Siggen_Setup *setup) {
  FILE *fp;
  char *cp, line[MAX_LINE], fname_unf[MAX_LINE];
  int i, j, k;
  double x, y, z;
  double fx, fy, fz;
  point pt;
  vector zero = {0,0,0};
  int lineno = 0;
  static vector ***efld = NULL;


  if (efld == NULL) {  // here we assume that setup->numx, setup->numy, setup->numz never
                       // change from their initial values, which is reasonable
    if ((efld = (vector ***) malloc(setup->numx*sizeof(*efld))) == NULL) {
      error("Malloc 1 failed in setup_efield\n");
      return 1;
    }
    for (i = 0; i < setup->numx; i++) {
      if ((efld[i] = (vector **) malloc(setup->numy*sizeof(*efld[i]))) == NULL) {
	error("Malloc 2 failed in setup_efield\n");
	return 1;
      }
      for (j = 0; j < setup->numy; j++) {
	if ((efld[i][j] = (vector *) malloc(setup->numz*sizeof(*efld[i][j]))) == NULL) {
	  error("Malloc 3 failed in setup_efield\n");
	  return 1;
	}
      }
    }
  }
  /* zero the table */
  for (i = 0; i < setup->numx; i++) {
    for (j = 0; j < setup->numy; j++) {
      for (k = 0; k < setup->numz; k++) {
	efld[i][j][k] = zero;
      }
    }
  }

  /* try to read from unformatted file */
  strncpy(fname_unf, setup->field_name, sizeof(fname_unf));
  if (!strstr(fname_unf, "unf")) strncat(fname_unf, "_unf", sizeof(fname_unf)-4);
  if ((fp = fopen(fname_unf, "r"))) {
    tell(NORMAL, "Reading field from unformatted file: %s\n", fname_unf);
    while (fgets(line, sizeof(line), fp) && line[0] == '#' &&
           !strstr(line, "start of unformatted data"))
      ;
    if (line[0] != '#') rewind(fp);
    for (i = 0; i < setup->numx; i++) {
      for (j = 0; j < setup->numy; j++) {
        if (fread(efld[i][j], sizeof(vector), setup->numz, fp) != setup->numz) {
          error("Error while reading %s\n", fname_unf);
          return -1;
        }
      }
    }
    fclose(fp);
    tell(NORMAL, "Done reading field\n");
    setup->efld = efld;
    return 0;
  }

  /* that didn't work; read from formatted file */
  if ((fp = fopen(setup->field_name, "r")) == NULL) {
    error("failed to open electric field table: %s\n", setup->field_name);
    return 1;
  }
  
  tell(NORMAL, "Reading electric field data from file: %s\n", setup->field_name);
  /* now read the table */
   while(fgets(line, MAX_LINE, fp) != NULL) {
    lineno++;
    for (cp = line; isspace(*cp) && *cp != '\0'; cp++);
    // discard comment lines
    if (*cp == '#' || !strlen(cp)) continue;

    if (sscanf(line, "%lf %lf %lf %lf %lf %lf", &x, &y, &z, 
	       &fx, &fy, &fz) != 6) {
      if (strstr(line, "No Solution") != NULL &&
	  sscanf(line, "%lf %lf %lf", &x, &y, &z) == 3) {
	fx = fy = fz = 0.0;
      } else {
	error("failed to read electric field data from line no %d\n"
	      "of file %s\n", lineno, setup->field_name);
	error("line is: %s\n", line);
	fclose(fp);
        return 1;
      }
    }
    i = lrint((x - setup->xmin)/setup->xtal_grid);
    j = lrint((y - setup->ymin)/setup->xtal_grid);
    k = lrint((z - setup->zmin)/setup->xtal_grid);
    if (i < 0 || i >= setup->numx || j < 0 || j >= setup->numy || k < 0 || k >= setup->numz)
      continue;
    pt.x = x; pt.y = y; pt.z = z;
    if (!in_crystal(setup, pt)) continue;
    efld[i][j][k].x = fx; efld[i][j][k].y = fy; efld[i][j][k].z = fz;
  }
  tell(NORMAL, "Done reading %d lines of electric field data\n", lineno);
  fclose(fp);

  setup->efld = efld;

  if (!strstr(setup->field_name, "unf") && (fp = fopen(fname_unf, "w"))) {
    tell(NORMAL, "Saving field to unformatted file: %s\n", fname_unf);
    for (i = 0; i < setup->numx; i++) {
      for (j = 0; j < setup->numy; j++) {
        if (fwrite(efld[i][j], sizeof(vector), setup->numz, fp) != setup->numz) {
          error("Error while writing %s\n", fname_unf);
        }
      }
    }
    fclose(fp);
  }
  return 0;
}

/*setup_wp
  read weighting potential values from files. returns 0 on success*/
static int setup_wp(MJD_Siggen_Setup *setup)
{
  FILE   *fp;
  char   *cp, line[MAX_LINE], fname[MAX_LINE], fname_unf[MAX_LINE];
  int    i, j, k;
  double x, y, z, wp;
  int    lineno;
  point  pt;
  static double ***wpot = NULL;


  if (wpot == NULL) {//assuming setup->numx, setup->numy, setup->numz never change as for setup_efld
    if ((wpot = malloc(setup->numx*sizeof(*wpot))) == NULL) {
      error("Malloc failed in setup_wp\n");
      return 1;
    }
    for (i = 0; i < setup->numx; i++) {
      if ((wpot[i] = malloc(setup->numy*sizeof(*wpot[i]))) == NULL) {
	error("Malloc failed in setup_wp\n");
	//NB: memory leak here.
	return 1;
      }
      for (j = 0; j < setup->numy; j++) {
	if ((wpot[i][j] = malloc(setup->numz*sizeof(*wpot[i][j]))) == NULL) {
	  error("Malloc failed in setup_wp\n");
	  //memory leak again
	  return 1;
	}
      }
    }
  }
  /*zero the wp table*/
  for (i = 0; i < setup->numx; i++) {
    for (j = 0; j < setup->numy; j++) {
      for (k = 0; k < setup->numz; k++) {
        wpot[i][j][k] = 0.0;
      }
    }
  }

  if (strstr(setup->wp_name, "unf")) {  //  file is unformatted
    double *wp2 = malloc(setup->numz*sizeof(double));
    tell(NORMAL, "Reading weighting potentials from file %s\n", setup->wp_name);
    if ((fp = fopen(setup->wp_name, "r")) == NULL) {
      error("failed to open file: %s\n", setup->wp_name);
      return -1;
    }

    while (fgets(line, sizeof(line), fp) && line[0] == '#' &&
           !strstr(line, "start of unformatted data"))
      ;
    if (line[0] != '#') rewind(fp);
    for (i = 0; i < setup->numx; i++) {
      for (j = 0; j < setup->numy; j++) {
        if (fread(wp2, sizeof(double), setup->numz, fp) != setup->numz) {
          error("Error while reading %s\n", fname);
          return -1;
        }
        for (k = 0; k < setup->numz; k++) wpot[i][j][k] = wp2[k];
      }
    }
    fclose(fp);
    tell(NORMAL, "Done reading weighting potentials: %d x %d x %d values\n",
         i, j, k);
    setup->wpot = wpot;
    return 0;
  }

  /* file is formatted */
  if ((fp = fopen(setup->wp_name, "r")) == NULL) {
    error("failed to open file: %s\n", setup->wp_name);
    return -1;
  }

  tell(NORMAL, "Reading weighting potential from file: %s\n", setup->wp_name);
  lineno = 0;
  while (fgets(line, MAX_LINE, fp) != NULL) {
    lineno++;
    /* discard comment lines */
    for (cp = line; isspace(*cp) && *cp != '\0'; cp++);
    if (*cp == '#' || !strlen(cp)) continue;
    if (sscanf(line, "%lf %lf %lf %lf\n",&x, &y, &z, &wp) != 4) { 
      if (strstr(line, "No Solution") != NULL
          && sscanf(line, "%lf %lf %lf", &x, &y, &z) == 3) {
        wp = 0.0;
      } else {
        error("failed to read weighting potential from line %d\n"
              "line: %s", lineno, line);
        fclose(fp);
        return 1;
      }
    }
    i = lrintf((x - setup->xmin)/setup->xtal_grid);
    j = lrintf((y - setup->ymin)/setup->xtal_grid);
    k = lrintf((z - setup->zmin)/setup->xtal_grid);
    if (i < 0 || i >= setup->numx ||
        j < 0 || j >= setup->numy ||
        k < 0 || k >= setup->numz) continue;
    pt.x = x; pt.y = y; pt.z = z;
    if (!in_crystal(setup, pt)) continue;
    wpot[i][j][k] = wp;
  }
  tell(NORMAL,"Read %d lines from file %s\n", lineno, setup->wp_name);
  fclose(fp);

  setup->wpot = wpot;

  strncpy(fname_unf, setup->wp_name, sizeof(fname_unf));
  strncat(fname_unf, "_unf", sizeof(fname_unf)-4);
  if (!strstr(setup->wp_name, "unf") && (fp = fopen(fname_unf, "w"))) {
    tell(NORMAL, "Saving WPs to unformatted file: %s\n", fname_unf);
    for (i = 0; i < setup->numx; i++) {
      for (j = 0; j < setup->numy; j++) {
        if (fwrite(wpot[i][j], sizeof(float), setup->numz, fp) != setup->numz) {
          error("Error while writing %s\n", fname_unf);
        }
      }
    }
    tell(NORMAL, "Done saving weighting potentials\n");
    fclose(fp);
  }
  return 0;
}

/* free malloc()'ed memory and do other cleanup*/
int fields_finalize(MJD_Siggen_Setup *setup) {
  int i, j;

  for (i = 0; i < setup->numx; i++) {
    for (j = 0; j < setup->numy; j++) {
      free(setup->efld[i][j]);
    }
    free(setup->efld[i]);
  }
  free(setup->efld);
  for (i = 0; i < setup->numx; i++) {
    for (j = 0; j < setup->numy; j++) {
      free(setup->wpot[i][j]);
    }
    free(setup->wpot[i]);
  }
  free(setup->wpot);
  free(setup->v_lookup);
 
  return 1;
}

void set_temp(MJD_Siggen_Setup *setup, float temp) {
  if (temp < MIN_TEMP || temp > MAX_TEMP) {
    printf("temperature out of range: %f\n", temp);
  } else {
    setup->xtal_temp = temp;
    printf("temperature set to %f\n", setup->xtal_temp);
  }
}
