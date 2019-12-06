/* program to calculate weighting potentials of coaxial Ge detectors
   by relaxation

   to run: ./fieldgen3d <config_file_name> -b bias_volts -w {0,1} -d {0,1} -p {0,1})
                        -z rho_spectrum_file_name

*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
float fminf(float x, float y);
#include <unistd.h>
#include <ctype.h>

#include "mjd_siggen.h"

#define e_over_E (1.413*8/3)



static int grid_init(MJD_Siggen_Setup *setup);
static int ev_calc(MJD_Siggen_Setup *setup, MJD_Siggen_Setup *old_setup);
static int wp_calc(MJD_Siggen_Setup *setup, MJD_Siggen_Setup *old_setup);
static int malloc_arrays(MJD_Siggen_Setup *setup);
static int write_ev(MJD_Siggen_Setup *setup, char *fname);
static int write_wp(MJD_Siggen_Setup *setup, char *fname);
int do_relaxation(MJD_Siggen_Setup *setup, int efld_calc);

int main(int argc, char *argv[]) {
  MJD_Siggen_Setup setup1, setup2, setup3, *setup = &setup1;
  float rho_z[256] = {0};
  int   i;
  FILE  *fp;


  if (argc%2 != 0 || !strncmp(argv[1], "-h", 2)) {
    printf("usage: %s <config_file_name>\n"
           "        -b bias_volts\n"
           "        -w {0,1}   (do_not/do write the field file)\n"
           "        -p {0,1}   (do_not/do write the WP file)\n"
           "        -f {0,1}   (do_not/do try to fix potential values from larger grid sizes)\n"
           "                    - note that this reduces accuracy of min field and depl voltage\n"
           "        -z rho_spectrum_file_name\n", argv[0]);
    return 0;
  }
  if (read_config(argv[1], &setup1)) return 1;
  strncpy(setup1.config_file_name, argv[1], 256);
  setup1.fix_adaptive = 0;

  for (i=2; i<argc-1; i+=2) {
    if (strstr(argv[i], "-b")) {
      setup1.xtal_HV = atof(argv[i+1]);     // bias volts
    } else if (strstr(argv[i], "-w")) {
      setup1.write_field = atoi(argv[i+1]); // write-out options
    } else if (strstr(argv[i], "-p")) {
      setup1.write_WP = atoi(argv[i+1]);    // weighting-potential options
    } else if (strstr(argv[i], "-f")) {
      setup1.fix_adaptive = atoi(argv[i+1]);      // adaptive-grid-potential fixing option
    } else if (strstr(argv[i], "-r")) {
      if (!(fp = fopen(argv[i+1], "r"))) { // impurity-profile-spectrum file name
        printf("\nERROR: cannot open impurity profile spectrum file %s\n\n", argv[i+1]);
        return 1;
      }
      fread(rho_z, 36, 1, fp);
      fread(rho_z, sizeof(rho_z), 1, fp);
      fclose(fp);
      printf(" z(mm)   rho\n");
      for (i=0; i < 200 && rho_z[i] != 0.0f; i++)  printf(" %3d  %7.3f\n", i, rho_z[i]);
    } else {
      printf("Possible options:\n"
	     "      -b bias_volts\n"
	     "      -w {0,1,2} (for WV options)\n"
	     "      -p {0,1}   (for WP options)\n"
             "      -f {0,1}   (do_not/do try to fix potential values from larger grid sizes)\n"
             "                  - note that this reduces accuracy of min field and depl voltage\n"
             "      -r rho_spectrum_file_name\n");
      return 1;
    }
    // FIXME - implement other option: -z
  }
  printf("write_WP = %d\n", setup1.write_WP);

  printf("number of iterations: %d\n", setup1.max_iterations);

  memcpy(&setup2, &setup1, sizeof(setup1));
  memcpy(&setup3, &setup1, sizeof(setup1));

  setup1.xtal_grid *= 4.0;
  setup2.xtal_grid *= 2.0;
  if (grid_init(&setup1) != 0 ||
      grid_init(&setup2) != 0 ||
      grid_init(&setup3) != 0) {
    error("failed to init field calculations\n");
    return 1;
  }
  if (setup->write_field) {
    setup1.write_field = 0; // no need to save intermediate calculations
    setup2.write_field = 0;
    if (setup3.xtal_grid > 0.4) {
    ev_calc(&setup2, NULL);
    } else {
      ev_calc(&setup1, NULL);
      ev_calc(&setup2, &setup1);
    }
    ev_calc(&setup3, &setup2);
  }
  if (setup->write_WP) {
    setup1.write_WP = 0; // no need to save intermediate calculations
    setup2.write_WP = 0;
    if (setup3.xtal_grid > 0.4) {
      wp_calc(&setup2, NULL);
    } else {
      wp_calc(&setup1, NULL);
      wp_calc(&setup2, &setup1);
    }
    wp_calc(&setup3, &setup2);
  }
 
  return 0;
}


static int grid_init(MJD_Siggen_Setup *setup) {

  setup->xmin = setup->ymin = - setup->xtal_radius;
  setup->xmax = setup->ymax = setup->xtal_radius;
  setup->zmin = 0; setup->zmax = setup->xtal_length;

  setup->numx = (int)rint((setup->xmax - setup->xmin)/setup->xtal_grid) + 1;
  setup->numy = (int)rint((setup->ymax - setup->ymin)/setup->xtal_grid) + 1;
  setup->numz = (int)rint((setup->zmax - setup->zmin)/setup->xtal_grid) + 1;
  printf("numx, numy, numz: %d %d %d\n", setup->numx, setup->numy, setup->numz);
  
  if (malloc_arrays(setup) < 0) return -1;

  if (geometry_init(setup) <0) {
    error("failed to init geometry\n");
    return 1;
  }

  if (setup->impurity_lamda != 0 || setup->rho_b != 0) {
    setup->rho_c = -setup->impurity_gradient/10 -
                    setup->rho_b / setup->xtal_length * (1 - exp(-setup->xtal_length / setup->impurity_lamda));
    printf("rho_c = %f\n", setup->rho_c);
  }

  printf("\n\n"
	 "          grid size: %5.2f\n"
	 "         xmin, xmax: %5.2f %5.2f\n"
	 "         ymin, ymax: %5.2f %5.2f\n"
	 "         zmin, zmax: %5.2f %5.2f\n"
	 "     HV, rho0, drho: %.0f %.4f %.4f\n"
	 "lamda, rho_b, rho_c: %.4f %.4f %.4f\n", 
	 setup->xtal_grid, setup->xmin, setup->xmax, setup->ymin, setup->ymax,
         setup->zmin, setup->zmax, setup->xtal_HV, setup->impurity_z0, setup->impurity_gradient,
         setup->impurity_lamda, setup->rho_b, setup->rho_c);
  printf("\n");

  return 0;
}


int do_relaxation(MJD_Siggen_Setup *setup, int efld_calc) {

  int   i, j, k, *ylo, *yhi;
  int   old, new, iter, undep;
  float sum_dif, max_dif, dif, save_dif;
  float mean, z, x, y;
  double or_fact;
  if (efld_calc) {
    if (setup->xtal_grid > 0.7) or_fact = 0.89;
    else if (setup->xtal_grid > 0.4) or_fact = 0.95;
    else if (setup->xtal_grid > 0.2) or_fact = 0.972;
    else or_fact = 0.975;
  } else {
    if (setup->xtal_grid > 0.7) or_fact = 0.90;
    else if (setup->xtal_grid > 0.4) or_fact = 0.94;
    else if (setup->xtal_grid > 0.2) or_fact = 0.962;
    else or_fact = 0.965;
  }

  ylo = malloc(setup->numx*sizeof(*ylo));
  yhi = malloc(setup->numx*sizeof(*yhi));

  old = 1; new = 0;
  for (i = 0; i < setup->numx; i++) {
    ylo[i] = setup->numy-1;
    yhi[i] = 1;
    for (j = 0; j < setup->numy; j++) {
      for (k = 0; k < setup->numz; k++) {
        setup->e[i][j][k] = 0;
	if (setup->point_type[i][j][k] == OUTSIDE) continue;
        if (ylo[i] > j) ylo[i] = j;
        if (yhi[i] < j) yhi[i] = j;
	if (setup->point_type[i][j][k] != INSIDE &&
            setup->point_type[i][j][k] != FIXED &&
            setup->point_type[i][j][k] != CONTACT_0 &&
            setup->point_type[i][j][k] != CONTACT_VB &&
            setup->point_type[i][j][k] != PASSIVE &&
            setup->point_type[i][j][k] != DITCH &&
	    (i > 0 && i + 1 < setup->numx) &&
            (j > 0 && j + 1 < setup->numy) &&
            (k > 0 && k < setup->numz)) {
          printf("error for point %d %d %d\n", i, j, k);
	}

	if ((i > 0             && setup->point_type[i-1][j][k] == OUTSIDE &&
	    i+1 < setup->numx  && setup->point_type[i+1][j][k] == OUTSIDE) ||
	    (j > 0             && setup->point_type[i][j-1][k] == OUTSIDE &&
             j+1 < setup->numy && setup->point_type[i][j+1][k] == OUTSIDE)) {
	  printf("error2 for point %d %d %d\n", i, j, k);
	}
      }
    }
  }

  /* calculate voxel impurity value as a function of z */
  for (k = 0; k < setup->numz; k++) {
    z = setup->zmin + k*setup->xtal_grid;
    if (setup->impurity_lamda == 0.0 || setup->rho_b == 0.0) {
      setup->impurity_z[k] = (setup->impurity_z0 + setup->impurity_gradient * z/10) *
        3.0 * setup->xtal_grid*setup->xtal_grid  * e_over_E/6.0;
    } else {
      setup->impurity_z[k] = (setup->impurity_z0 -
                                    setup->rho_b*(1-exp(-z/setup->impurity_lamda)) -
                                    setup->rho_c*z) *
        3.0 * setup->xtal_grid*setup->xtal_grid * e_over_E/6.0;
    }
  }
  /* add surface charge to passivated surface */
  setup->impurity_z[setup->numz-1] +=
    setup->impurity_surface * 3.0 * setup->xtal_grid * e_over_E/6.0;

  printf("starting relaxation...\n");
  for (iter = 0; iter < setup->max_iterations; iter++) {
    old = 1 - old;
    new = 1 - new;
    sum_dif = 0.0;
    max_dif = 0.0;

    /* do relaxation iteration */
    for (i = 1; i < setup->numx-1; i++) {
      x = setup->xmin + i*setup->xtal_grid;
      for (j = ylo[i]; j <= yhi[i]; j++) {
	y = setup->ymin + j*setup->xtal_grid;
	for (k = 0; k < setup->numz-1; k++) {
          if (setup->point_type[i][j][k] == OUTSIDE ||
              setup->point_type[i][j][k] == FIXED ||
              setup->point_type[i][j][k] == CONTACT_0 ||
              setup->point_type[i][j][k] == CONTACT_VB ) continue;
	  z = setup->zmin + k*setup->xtal_grid;
          // save step difference from previous iteration
          save_dif = setup->v[old][i][j][k] - setup->v[new][i][j][k];
          setup->e[i][j][k] = 0.5*setup->e[i][j][k] + 0.5*save_dif;  // for possible exponential over-relax value
          if (iter < 2) save_dif = setup->e[i][j][k] = 0; 
	  if (setup->point_type[i][j][k] == INSIDE) {                // main bulk volume
            mean = (setup->v[old][i-1][j][k] +
                    setup->v[old][i+1][j][k] +
                    setup->v[old][i][j-1][k] +
                    setup->v[old][i][j+1][k] +
                    setup->v[old][i][j][k-1] +
                    setup->v[old][i][j][k+1]);
	  } else if (setup->point_type[i][j][k] == PASSIVE) {        // passivated surface
            mean = (setup->v[old][i-1][j][k] +
                    setup->v[old][i+1][j][k] +
                    setup->v[old][i][j-1][k] +
                    setup->v[old][i][j+1][k] +
                    2.0*setup->v[old][i][j][k+1]);
	  } else if (setup->point_type[i][j][k] == DITCH) {          // ditch
            mean = (setup->v[old][i-1][j][k] * setup->epsilon[i-1][j][k] +
                    setup->v[old][i+1][j][k] * setup->epsilon[i+1][j][k] +
                    setup->v[old][i][j-1][k] * setup->epsilon[i][j-1][k] +
                    setup->v[old][i][j+1][k] * setup->epsilon[i][j+1][k] +
                    setup->v[old][i][j][k+1] * setup->epsilon[i][j][k+1]);
            if (k > 0) {
              mean += setup->v[old][i][j][k-1] * setup->epsilon[i][j][k-1];
              mean *= 6.0/(setup->epsilon[i-1][j][k] + setup->epsilon[i+1][j][k] +
                           setup->epsilon[i][j-1][k] + setup->epsilon[i][j+1][k] +
                           setup->epsilon[i][j][k+1] + setup->epsilon[i][j][k-1]);                           
            } else {
              mean += setup->v[old][i][j][k+1] * setup->epsilon[i][j][k+1];
              mean *= 6.0/(setup->epsilon[i-1][j][k] + setup->epsilon[i+1][j][k] +
                           setup->epsilon[i][j-1][k] + setup->epsilon[i][j+1][k] +
                           setup->epsilon[i][j][k+1] + setup->epsilon[i][j][k+1]);                           
            }
	  }

          setup->v[new][i][j][k] = mean/6.0;
          if (efld_calc && setup->point_type[i][j][k] != DITCH)
            setup->v[new][i][j][k] += setup->impurity_z[k];
	  dif = setup->v[old][i][j][k] - setup->v[new][i][j][k];
 
          //setup->v[new][i][j][k] += or_fact * save_dif;   // do over-relaxation
          setup->v[new][i][j][k] += or_fact * setup->e[i][j][k];;   // do over-relaxation

	  if (dif < 0.0) dif = -dif;
          /*
            float dif2 = setup->v[old][i][j][k] - setup->v[new][i][j][k];
            if (dif2 < 0.0) dif2 = -dif2;
            if (dif > dif2) dif = dif2;
          */
	  sum_dif += dif;
	  if (max_dif < dif) max_dif = dif;
	}
      }
    }
    sum_dif /= setup->numx * setup->numy * setup->numz;
    if (iter %100 == 0 ||
        (setup->xtal_grid < 0.4 && iter %50 == 0)) {
      printf("%5d %.3e %.3e\n", iter, max_dif, sum_dif);
      fflush(stdout);
    }
    if (max_dif < 1.0e-8) break;

    if (efld_calc) {
      //if (max_dif/setup->xtal_HV < 2.0e-3/5000.0 &&
      //    sum_dif/setup->xtal_HV < 3.0e-4/5000.0) break;
      if (max_dif < 1e-3 &&
          sum_dif < 6e-5) break;
    } else {
      //if (max_dif < 5.0e-7 &&
      //    sum_dif < 1.0e-7) break;
      if (max_dif < 3.0e-7 &&
          sum_dif < 2.0e-9) break;
    }
  }
  printf("%5d %.3e %.3e\n", 
	 iter, max_dif, sum_dif);

  undep = 0;
  for (i = 0; i < setup->numx && !undep; i++) {
    for (j = 0; j < setup->numy && !undep; j++) {
      for (k = 0; k < setup->numz && !undep; k++) {
	if (setup->v[new][i][j][k] > setup->xtal_HV ||
            setup->v[new][i][j][k] < 0) {
	  printf("detector is undepleted!\n");
	  undep = 1;
	}
      }
    }
  }
  printf("...  relaxation done\n");
  free(ylo); free(yhi);

  return 0;
}

static int v_interpolate(MJD_Siggen_Setup *setup, MJD_Siggen_Setup *old_setup) {
  int   i, j, k, i2, j2, k2, n1=0, n2=0;
  float fx, fy, fz, thresh;
  double ***ov = old_setup->v[1];

  for (i = 2; i < setup->numx-2; i++) {
    i2 = i/2;
    fx = 1.0 - 0.5 * (float) (i%2);
    for (j = 2; j < setup->numy-2; j++) {
      j2 = j/2;
      fy = 1.0 - 0.5 * (float) (j%2);
      for (k = 2; k < setup->numz-2; k++) {
        k2 = k/2;
        if (setup->point_type[i][j][k] >= INSIDE &&
            old_setup->point_type[i2][j2][k2] >= INSIDE) {
          fz = 1.0 - 0.5 * (float) (k%2);
          setup->v[0][i][j][k] = setup->v[1][i][j][k] =
            fx       * fy       * fz       * ov[i2  ][j2  ][k2  ] +
            (1.0-fx) * fy       * fz       * ov[i2+1][j2  ][k2  ] +
            fx       * (1.0-fy) * fz       * ov[i2  ][j2+1][k2  ] +
            (1.0-fx) * (1.0-fy) * fz       * ov[i2+1][j2+1][k2  ] +
            fx       * fy       * (1.0-fz) * ov[i2  ][j2  ][k2+1] +
            (1.0-fx) * fy       * (1.0-fz) * ov[i2+1][j2  ][k2+1] +
            fx       * (1.0-fy) * (1.0-fz) * ov[i2  ][j2+1][k2+1] +
            (1.0-fx) * (1.0-fy) * (1.0-fz) * ov[i2+1][j2+1][k2+1];
          n1++;
          thresh = fabs(ov[i2][j2][k2]) * 0.005;
          if (setup->fix_adaptive &&
              /*
              old_setup->point_type[i/2-1][j/2-1][k/2-1] >= INSIDE &&
              old_setup->point_type[i/2+1][j/2-1][k/2-1] >= INSIDE &&
              old_setup->point_type[i/2-1][j/2+1][k/2-1] >= INSIDE &&
              old_setup->point_type[i/2+1][j/2+1][k/2-1] >= INSIDE &&
              old_setup->point_type[i/2-1][j/2-1][k/2+1] >= INSIDE &&
              old_setup->point_type[i/2+1][j/2-1][k/2+1] >= INSIDE &&
              old_setup->point_type[i/2-1][j/2+1][k/2+1] >= INSIDE &&
              old_setup->point_type[i/2+1][j/2+1][k/2+1] >= INSIDE &&
              */
              fabs(ov[i2-1][j2][k2] + ov[i2+1][j2][k2] - 2.0*ov[i2][j2][k2]) < thresh &&
              fabs(ov[i2][j2-1][k2] + ov[i2][j2+1][k2] - 2.0*ov[i2][j2][k2]) < thresh &&
              fabs(ov[i2][j2][k2-1] + ov[i2][j2][k2+1] - 2.0*ov[i2][j2][k2]) < thresh) {
            setup->point_type[i][j][k] = FIXED;
            n2++;
          }
        }
      }
    }
  }
  printf("n1 = %d n2 = %d (%.4f)\n", n1, n2, (float) n2 / (float) n1);
  return 0;
}

static int ev_calc(MJD_Siggen_Setup *setup, MJD_Siggen_Setup *old_setup) {

  if (!old_setup) printf("\n\n ---- starting EV calculation --- \n");
  init_ev_calc(setup);
  if (old_setup) v_interpolate(setup, old_setup);
  do_relaxation(setup, 1);
  if (setup->write_field) {
    printf("saving E, V in file %s\n", setup->field_name);
    write_ev(setup, setup->field_name);
  }
  return 0;
}

static int wp_calc(MJD_Siggen_Setup *setup, MJD_Siggen_Setup *old_setup) {

  if (!old_setup) printf("\n\n ---- starting WP calculation --- \n");
  init_wp_calc(setup);
  if (old_setup) v_interpolate(setup, old_setup);
  do_relaxation(setup, 0);
  if (setup->write_WP) {
    printf("saving WP in file %s\n", setup->wp_name);
    write_wp(setup, setup->wp_name);
  }
  return 0;
}

int report_config(FILE *fp_out, char *config_file_name) {

  char  *c, line[256];
  FILE  *file;

  fprintf(fp_out, "# Config file: %s\n", config_file_name);
  if (!(file = fopen(config_file_name, "r"))) return 1;

  while (fgets(line, sizeof(line), file)) {
    if (strlen(line) < 3 || *line == ' ' || *line == '\t' || *line == '#') continue;
    if ((c = strchr(line, '#')) || (c = strchr(line, '\n'))) *c = '\0';
    fprintf(fp_out, "# %s\n", line);
  }
  fclose(file);
  return 0;
}

int write_ev(MJD_Siggen_Setup *setup, char *fname) {
  int i,j,k;
  FILE *fp = NULL, *fp2 = NULL;
  struct efld{float x; float y; float z;} ***e;


  if ((fp = fopen(fname, "w")) == NULL) {
    error("Failed to open output file: %s\n", fname);
    return -1;
  }
  if ((fp2 = fopen("fields/ev2d.dat", "w")) == NULL) {
    error("Failed to open output file: %s\n", fname);
    return -1;
  }
  /* copy configuration parameters to output file */
  report_config(fp, setup->config_file_name);
  fprintf(fp, "#\n# HV bias in fieldgen: %.1f V\n", setup->xtal_HV);

  if ((e = malloc(setup->numx*sizeof(*e))) == NULL) {
    error("malloc failed\n");
    fclose(fp);
    return -1;
  }
  for (i = 0; i < setup->numx; i++) {
    if ((e[i] = malloc(setup->numy*sizeof(*e[i]))) == NULL) {
      error("malloc failed\n");
      fclose(fp);
      return -1;
    }
    for (j = 0; j < setup->numy; j++) {
      if ((e[i][j] = calloc(setup->numz, sizeof(*e[i][j]))) == NULL) {
	error("malloc failed\n");
	fclose(fp);
	return -1;
      }
    }
  }
  for (i = 0; i < setup->numx; i++) {
    for (j = 0; j < setup->numy; j++) {
      for (k = 0; k < setup->numz; k++) {
	if (setup->point_type[i][j][k] == OUTSIDE) continue;
	if (setup->point_type[i][j][k] == CONTACT_0) continue;
	if (setup->point_type[i][j][k] == CONTACT_VB) continue;
	if (i == 0 || (setup->point_type[i-1][j][k] == OUTSIDE))
	  e[i][j][k].x = -10.0*(setup->v[0][i+1][j][k] - setup->v[0][i][j][k])/setup->xtal_grid;
	else if (i == setup->numx -1 || (setup->point_type[i+1][j][k] == OUTSIDE))
	  e[i][j][k].x = -10.0*(setup->v[0][i][j][k] - setup->v[0][i-1][j][k])/setup->xtal_grid;
	else
	  e[i][j][k].x = -10.0*(setup->v[0][i+1][j][k] - setup->v[0][i-1][j][k])/setup->xtal_grid/2;
	if (j == 0 || (setup->point_type[i][j-1][k] == OUTSIDE))
	  e[i][j][k].y = -10.0*(setup->v[0][i][j+1][k] - setup->v[0][i][j][k])/setup->xtal_grid;
	else if (j == setup->numy -1 || setup->point_type[i][j+1][k] == OUTSIDE)
	  e[i][j][k].y = -10.0*(setup->v[0][i][j][k] - setup->v[0][i][j-1][k])/setup->xtal_grid;
	else
	e[i][j][k].y = -10.0*(setup->v[0][i][j+1][k] - setup->v[0][i][j-1][k])/setup->xtal_grid/2;
	if (k == 0 || (setup->point_type[i][j][k-1] == OUTSIDE))
	  e[i][j][k].z = -10.0*(setup->v[0][i][j][k+1] - setup->v[0][i][j][k])/setup->xtal_grid;
	else if (k == setup->numz -1 || (setup->point_type[i][j][k+1] == OUTSIDE))
	  e[i][j][k].z = -10.0*(setup->v[0][i][j][k] - setup->v[0][i][j][k-1])/setup->xtal_grid;
	else
	  e[i][j][k].z = -10.0*(setup->v[0][i][j][k+1] - setup->v[0][i][j][k-1])/setup->xtal_grid/2;
      }
    }
  }
  if (strstr(fname, "unf")) {
    fprintf(fp, "#\n## start of unformatted data\n");
    for (i = 0; i < setup->numx; i++) {
      for (j = 0; j < setup->numy; j++) {
        if (fwrite(e[i][j], sizeof(vector), setup->numz, fp) != setup->numz) {
          error("Error while writing %s\n", fname);
        }
      }
    }
  } else {
    fprintf(fp, "#\n## x[mm]  y[mm]  z[mm]  Ex[V/cm]  Ey[V/cm]  Ez[V/cm]  |E|[V/cm] V[V]\n");
#define SQ(x) ((x)*(x))
    for (i = 0; i < setup->numx; i++) {
      for (j = 0; j < setup->numy; j++) {
        for (k = 0; k < setup->numz; k++) {
          fprintf(fp, "%6.2f %6.2f %6.2f\t %10e %10e %10e %10e  %10e\n",
                  (setup->xmin+i*setup->xtal_grid), 
                  (setup->ymin+j*setup->xtal_grid),
                  (setup->zmin+k*setup->xtal_grid), 		
                  e[i][j][k].x, e[i][j][k].y, e[i][j][k].z,
                  sqrt(SQ(e[i][j][k].x) + SQ(e[i][j][k].y) + SQ(e[i][j][k].z)),
                  setup->v[0][i][j][k]);
        }
      }
    }
  }

  fprintf(fp2, "#\n## x[mm]  y[mm]  z[mm]  Ex[V/cm]  Ey[V/cm]  Ez[V/cm]  |E|[V/cm] V[V]\n");
  i = setup->numx/2;
  for (j = 0; j < setup->numy; j++) {
    for (k = 0; k < setup->numz; k++) {
      fprintf(fp2, "%6.2f %6.2f %6.2f\t %10e %10e %10e %10e  %10e\n",
              (setup->xmin+i*setup->xtal_grid), 
              (setup->ymin+j*setup->xtal_grid),
              (setup->zmin+k*setup->xtal_grid), 		
              e[i][j][k].x, e[i][j][k].y, e[i][j][k].z,
              sqrt(SQ(e[i][j][k].x) + SQ(e[i][j][k].y) + SQ(e[i][j][k].z)),
              setup->v[0][i][j][k]);
    }
  }
  for (i = 0; i < setup->numx; i++) {
    for (j = 0; j < setup->numy; j++) { 
      free(e[i][j]);
    }
    free(e[i]);
  }
  free(e);
  
  fclose(fp);
  fclose(fp2);
  return 0;
 }

 int write_wp(MJD_Siggen_Setup *setup, char *fname) {
  int i,j,k;
  FILE *fp = NULL, *fp2 = NULL;

  if ((fp = fopen(fname, "w")) == NULL) {
    error("Failed to open output file: %s\n", fname);
    return -1;
  }
  if ((fp2 = fopen("fields/wp2d.dat", "w")) == NULL) {
    error("Failed to open output file: %s\n", fname);
    return -1;
  }
  /* copy configuration parameters to output file */
  report_config(fp, setup->config_file_name);

  if (strstr(fname, "unf")) {
    fprintf(fp, "#\n## start of unformatted data\n");
    for (i = 0; i < setup->numx; i++) {
      for (j = 0; j < setup->numy; j++) {
        if (fwrite(setup->v[0][i][j], sizeof(double), setup->numz, fp) != setup->numz) {
          error("Error while writing %s\n", fname);
        }
      }
    }
  } else {
    fprintf(fp, "#\n## x[mm]  y[mm]  z[mm]   WP\n");
    for (i = 0; i < setup->numx; i++) {
      for (j = 0; j < setup->numy; j++) {
        for (k = 0; k < setup->numz; k++) {
          fprintf(fp, "%6.2f %6.2f %6.2f\t %10e\n",
                  (setup->xmin + i * setup->xtal_grid),
                  (setup->ymin + j * setup->xtal_grid),
                  (setup->zmin + k * setup->xtal_grid), 
                  setup->v[0][i][j][k]);
        }
      }
    }
  }

  fprintf(fp2, "#\n## x[mm]  y[mm]  z[mm]   WP\n");
  i = setup->numx/2;
  for (j = 0; j < setup->numy; j++) {
    for (k = 0; k < setup->numz; k++) {
      fprintf(fp2, "%6.2f %6.2f %6.2f\t %10e\n",
              (setup->xmin + i * setup->xtal_grid),
              (setup->ymin + j * setup->xtal_grid),
              (setup->zmin + k * setup->xtal_grid), 
              setup->v[0][i][j][k]);
    }
  }
  
  fclose(fp);
  fclose(fp2);
  return 0;
}


static int malloc_arrays(MJD_Siggen_Setup *setup) {
  int i, j, k, kk;
  kk = setup->ditch_depth/setup->xtal_grid + 2;

  if ((setup->impurity_z = malloc(setup->numz * sizeof(double))) == NULL ||
      (setup->e          = malloc(setup->numx * sizeof(setup->e[0]))) == NULL ||
      (setup->epsilon    = malloc(setup->numx * sizeof(setup->epsilon[0]))) == NULL ||
      (setup->v[0]       = malloc(setup->numx * sizeof(setup->v[0][0]))) == NULL ||
      (setup->v[1]       = malloc(setup->numx * sizeof(setup->v[0][0]))) == NULL ||
      (setup->point_type = malloc(setup->numx * sizeof(setup->point_type))) == NULL) {
    error("malloc failed\n");
    return -1;
  }
  for (i = 0; i < setup->numx; i++) {
    if ((setup->e[i]          = malloc(setup->numy * sizeof(*setup->e[i]))) == NULL ||
        (setup->epsilon[i]    = malloc(setup->numy * sizeof(*setup->epsilon[i]))) == NULL ||
        (setup->v[0][i]       = malloc(setup->numy * sizeof(*setup->v[0][i]))) == NULL ||
	(setup->v[1][i]       = malloc(setup->numy * sizeof(*setup->v[1][i]))) == NULL ||
	(setup->point_type[i] = malloc(setup->numy * sizeof(setup->point_type[i]))) == NULL) {
      error("malloc failed\n");
      return -1;
    }
    for (j = 0; j < setup->numy; j++) {
      if ((setup->e[i][j]          = malloc(setup->numz * sizeof(*setup->e[i][j]))) == NULL ||
          (setup->epsilon[i][j]    = malloc(kk * sizeof(*setup->epsilon[i][j]))) == NULL ||
          (setup->v[0][i][j]       = malloc(setup->numz * sizeof(*setup->v[0][i][j]))) == NULL ||
	  (setup->v[1][i][j]       = malloc(setup->numz * sizeof(*setup->v[1][i][j]))) == NULL ||
	  (setup->point_type[i][j] = calloc(setup->numz,  sizeof(setup->point_type[i][j]))) == NULL) {
	error("malloc failed\n");
	return -1;
      }
      for (k = 0; k < kk; k++) setup->epsilon[i][j][k] = 16.0;
    }
  }

  return 0;
}
