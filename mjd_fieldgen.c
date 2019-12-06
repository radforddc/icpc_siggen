/* program to calculate electric fields and weighting potentials
              of PPC BEGe, and inverted-coax-PC Ge detectors by relaxation
   author:           D.C. Radford
   first written:         Nov 2007
   MJD modified version:  Oct 2014, June 2016
      - uses the same (single) config file as modified MJD siggen
      - added intelligent coarse grid / refinement of grid
      - added interpolation of RC and LC positions on the grid
      - added optional bulletization of point contact
   This modified version: Dec 2016
      - added hole, inner taper length, outer taper length, taper angle
   Nov  2017: added top bulletization

   May-July 2019
      - added dead layer / Li thickness
      - added over-relaxation factor; gained factor of ~ 40 speedup, depending on size/grid
      - added code to estimate full-depletion voltage in the case of bubble depletion (pinch-off)

   Nov 2019
      - major rewrite and re-arragement of the code
      - broke out functions grid_init, ev_calc, wp_calc, write_ev, write_wp, do_relax, interpolate
      - got an additional factor of ~3 speedup
      - more modular and readable
      - todo now: 
          - add back in interpolation of point contact
          - execute WD=1 option

   TO DO:
      - on coarse grids, interpolate the position of L, R, LT, and the ditch
            (as is done now already for RC and LC)
*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "mjd_siggen.h"

#define MAX_ITS 50000     // default max number of iterations for relaxation
#define MAX_ITS_FACTOR 2  // factor by which max iterations is reduced as grid is refined

static int report_config(FILE *fp_out, char *config_file_name);
static int grid_init(MJD_Siggen_Setup *setup);
static int ev_calc(MJD_Siggen_Setup *setup, MJD_Siggen_Setup *old_setup);
static int wp_calc(MJD_Siggen_Setup *setup, MJD_Siggen_Setup *old_setup);
static int write_ev(MJD_Siggen_Setup *setup);
static int write_wp(MJD_Siggen_Setup *setup);
static int do_relax(MJD_Siggen_Setup *setup, int ev_calc);
static int ev_relax_undep(MJD_Siggen_Setup *setup);
static int wp_relax_undep(MJD_Siggen_Setup *setup);
static int interpolate(MJD_Siggen_Setup *setup, MJD_Siggen_Setup *old_setup);

/* -------------------------------------- main ------------------- */
int main(int argc, char **argv)
{

  MJD_Siggen_Setup setup, setup1, setup2;

  float BV;      // bias voltage
  int   WV = 0;  // 0: do not write the V and E values to ppc_ev.dat
                 // 1: write the V and E values to ppc_ev.dat
                 // 2: write the V and E values for both +r, -r (for gnuplot, NOT for siggen)
  int   WP = 0;  // 0: do not calculate the weighting potential
                 // 1: calculate the WP and write the values to ppc_wp.dat
  int   WD = 0;  // 0: do not write out depletion surface
                 // 1: write out depletion surface to depl_<HV>.dat

  int   i, j;
  FILE  *fp;


  if (argc < 2 || argc%2 != 0 || read_config(argv[1], &setup)) {
    printf("Usage: %s <config_file_name> [options]\n"
           "   Possible options:\n"
	   "      -b bias_volts\n"
	   "      -w {0,1}  do_not/do write the field file)\n"
	   "      -d {0,1}  do_not/do write the depletion surface)\n"
	   "      -p {0,1}  do_not/do write the WP file)\n"
           "      -r rho_spectrum_file_name\n", argv[0]);
    return 1;
  }
  strncpy(setup.config_file_name, argv[1], sizeof(setup.config_file_name));

  if (setup.xtal_grid < 0.001) setup.xtal_grid = 0.5;
  BV = setup.xtal_HV;
  WV = setup.write_field;
  WP = setup.write_WP;
  setup.rho_z_spe[0] = 0;

  for (i=2; i<argc-1; i++) {
    if (strstr(argv[i], "-b")) {
      BV = setup.xtal_HV = atof(argv[++i]);   // bias volts
    } else if (strstr(argv[i], "-w")) {
      WV = atoi(argv[++i]);   // write-out options
    } else if (strstr(argv[i], "-d")) {
      WD = atoi(argv[++i]);   // write-out options
    } else if (strstr(argv[i], "-p")) {
      WP = atoi(argv[++i]);   // weighting-potential options
    } else if (strstr(argv[i], "-r")) {
      if (!(fp = fopen(argv[++i], "r"))) {   // impurity-profile-spectrum file name
        printf("\nERROR: cannot open impurity profile spectrum file %s\n\n", argv[i+1]);
        return 1;
      }
      fread(setup.rho_z_spe, 36, 1, fp);
      for (j=0; j<1024; j++) setup.rho_z_spe[i] = 0;
      fread(setup.rho_z_spe, sizeof(setup.rho_z_spe), 1, fp);
      fclose(fp);
      printf(" z(mm)   rho\n");
      for (j=0; j < 200 && setup.rho_z_spe[j] != 0.0f; j++)
        printf(" %3d  %7.3f\n", j, setup.rho_z_spe[j]);
    } else {
      printf("Possible options:\n"
	     "      -b bias_volts\n"
	     "      -w {0,1,2} (for WV options)\n"
	     "      -p {0,1}   (for WP options)\n"
             "      -r rho_spectrum_file_name\n");
      return 1;
    }
  }

  if (setup.xtal_length/setup.xtal_grid * setup.xtal_radius/setup.xtal_grid > 2500*2500) {
    printf("Error: Crystal size divided by grid size is too large!\n");
    return 1;
  }
  if (WV < 0 || WV > 2) WV = 0;

  /* -------------- give details of detector geometry */
  if (setup.verbosity >= CHATTY) {
    printf("\n\n"
           "      Crystal: Radius x Length: %.1f x %.1f mm\n",
	   setup.xtal_radius, setup.xtal_length);
    if (setup.hole_length > 0) {
      if (setup.inner_taper_length > 0)
        printf("    Core hole: Radius x length: %.1f x %.1f mm,"
               " taper %.1f x %.1f mm (%2.f degrees)\n",
               setup.hole_radius, setup.hole_length,
               setup.inner_taper_width, setup.inner_taper_length, setup.taper_angle);
      else
        printf("    Core hole: Radius x length: %.1f x %.1f mm\n",
               setup.hole_radius, setup.hole_length);
    }
    printf("Point contact: Radius x length: %.1f x %.1f mm\n",
           setup.pc_radius, setup.pc_length);
    if (setup.ditch_depth > 0) {
      printf("  Wrap-around: Radius x ditch x gap:  %.1f x %.1f x %.1f mm\n",
             setup.wrap_around_radius, setup.ditch_depth, setup.ditch_thickness);
    }
    printf("         Bias: %.0f V\n", BV);
  }
    
  if ((BV < 0 && setup.impurity_z0 < 0) || (BV > 0 && setup.impurity_z0 > 0)) {
    printf("ERROR: Expect bias and impurity to be opposite sign!\n");
    return 1;
  } 
  if (setup.impurity_z0 > 0) {
    // swap polarity for n-type material; this lets me assume all voltages are positive
    BV = -BV;
    setup.xtal_HV *= -1.0;
    setup.impurity_z0         *= -1.0;
    setup.impurity_gradient   *= -1.0;
    setup.impurity_quadratic  *= -1.0;
    setup.impurity_surface    *= -1.0;
    setup.impurity_radial_add *= -1.0;
  }
  /* use an adaptive grid; start out coarse and then refine the grid */
  memcpy(&setup1, &setup, sizeof(setup));
  memcpy(&setup2, &setup, sizeof(setup));
  setup1.xtal_grid *= 9.0;
  setup2.xtal_grid *= 3.0;
  if (grid_init(&setup1) != 0 ||
      grid_init(&setup2) != 0 ||
      grid_init(&setup)  != 0) {
    printf("failed to init field calculations\n");
    return 1;
  }

  /* -------------- calculate electric potential/field */
  if (setup.write_field) {
    setup1.write_field = 0; // no need to save intermediate calculations
    setup2.write_field = 0;
    if (setup.xtal_grid > 0.4) {
      ev_calc(&setup2, NULL);
    } else {
      ev_calc(&setup1, NULL);
      ev_calc(&setup2, &setup1);
    }
    ev_calc(&setup, &setup2);
  }

  /* -------------- calculate weighting potential */
  if (setup.write_WP) {
    setup1.write_WP = 0; // no need to save intermediate calculations
    setup2.write_WP = 0;
    if (setup.xtal_grid > 0.4) {
      wp_calc(&setup2, NULL);
    } else {
      wp_calc(&setup1, NULL);
      wp_calc(&setup2, &setup1);
    }
    wp_calc(&setup, &setup2);
  }

  /* -------------- calculate capacitance
     1/2 * epsilon * integral(E^2) = 1/2 * C * V^2
     so    C = epsilon * integral(E^2) / V^2    V = 1 volt
  */
  double esum, esum2, pi=3.14159, Epsilon=(8.85*16.0/1000.0);  // permittivity of Ge in pF/mm
  float  E_r, E_z;
  float  grid = setup.xtal_grid;
  int    r, z, test;
  int    L  = lrint(setup.xtal_length/grid)+1;
  int    R  = lrint(setup.xtal_radius/grid)+1;
  int    LC = lrint(setup.pc_length/grid)+1;
  int    RC = lrint(setup.pc_radius/grid)+1;

  esum = esum2 = test = 0;
  for (z=1; z<L; z++) {
    for (r=2; r<R; r++) {
      E_r = setup.eps_dr[z][r]/16.0 * (setup.v[1][z][r] - setup.v[1][z][r+1])/(0.1*grid);
      E_z = setup.eps_dz[z][r]/16.0 * (setup.v[1][z][r] - setup.v[1][z+1][r])/(0.1*grid);
      esum += (E_r*E_r + E_z*E_z) * (double) (r-1);
      if ((r == RC   && z <= LC)   || (r <= RC   && z == LC)   ||
	  (r == RC+1 && z <= LC+1) || (r <= RC+1 && z == LC+1)) { // average over two different surfaces
	if (setup.point_type[z+1][r+1] == PC) test = 1;
	esum2 += 0.5 * sqrt(E_r*E_r + E_z*E_z) * (double) (r-1);  // 0.5 since averaging over 2 surfaces
      }
    }
  }
  esum  *= 2.0 * pi * 0.01 * Epsilon * pow(grid, 3.0);
  // Epsilon is in pF/mm
  // 0.01 converts (V/cm)^2 to (V/mm)^2, pow() converts to grid^3 to mm3
  esum2 *= 2.0 * pi * 0.1 * Epsilon * pow(grid, 2.0);
  // 0.1 converts (V/cm) to (V/mm),  grid^2 to  mm2
  printf("  >>  Calculated capacitance at %.0f V: %.3lf pF\n", BV, esum);
  if (!test)
    printf("  >>  Alternative calculation of capacitance: %.3lf pF\n", esum2);

  /* -------------- estimate depletion voltage */
  double min = BV, min2 = BV, dV, dW, testv;
  int    vminr=0, vminz=0;
  int    dz[4] = {1, -1, 0, 0}, dr[4] = {0, 0, 1, -1};
  if (setup.fully_depleted) {
    // find minimum potential
    for (z=1; z<LC+2; z++) {
      for (r=1; r<RC+2; r++) {
        if (setup.vsave[z][r] > 0 &&
            min > setup.vsave[z][r] / (1.0 - setup.v[1][z][r])) {
          min = setup.vsave[z][r] / (1.0 - setup.v[1][z][r]);
        }
      }
    }

    /* check for bubble depletion / pinch-off by seeing how much the bias
       must be reduced for any pixel to be in a local potential minimum  */
    for (z=LC+2; z<L-2; z++) {
      for (r=1; r<R-2; r++) {
        if (setup.point_type[z][r] == INSIDE && setup.v[1][z][r] > 0.0001) {
          testv = -1;
          for (i=0; i<4; i++) {
            if (r==1 && i==2) break;  // do not check dr for r=1 (=0.0)
            dV = setup.vsave[z+dz[i]][r+dr[i]]  - setup.vsave[z][r];  // potential
            dW = setup.v[1][z+dz[i]][r+dr[i]]   - setup.v[1][z][r];   // WP
            if (dW*grid > 0.00001 && dV < 0 && testv < -dV/dW) testv = -dV/dW;
          }
          if (testv >= 0 && min2 > testv) {
            min2 = testv;
            vminr = r; vminz = z;
          }
        }
      }
    }
    if (min2 < min) {
      printf("Estimated pinch-off voltage = %.0f V\n", BV - min);
      printf(" min2 = %.1f at (r,z) = (%.1f, %.1f), so\n",
             min2, (vminr-1)*grid, (vminz-1)*grid);
      printf("   Full depletion (max pinch-off voltage) = %.0f\n", BV - min2);
    } else {
      printf("Estimated depletion voltage = %.0f V\n", BV - min);
    }
    printf("Minimum bulk field = %.1f V/cm at (r,z) = (%.1f, %.1f) mm\n\n",
           setup.Emin, setup.rmin, setup.zmin);
  }
 
  return 0;
} /* main */

/* -------------------------------------- ev_calc ------------------- */
int ev_calc(MJD_Siggen_Setup *setup, MJD_Siggen_Setup *old_setup) {
  int    i, j;
  float  grid = setup->xtal_grid;
  int    L  = lrint(setup->xtal_length/grid)+3;
  int    R  = lrint(setup->xtal_radius/grid)+3;

  if (!old_setup) {
    printf("\n\n ---- starting EV calculation --- \n");
    for (i = 1; i < L; i++) {
      for (j = 1; j < R; j++) {
        setup->v[0][i][j] = setup->v[1][i][j] = setup->xtal_HV/2.0;
      }
    }
  }
  if (old_setup) interpolate(setup, old_setup);
  setup->fully_depleted = 1;
  setup->bubble_volts = 0;

  /* set boundary voltages */
  for (i = 1; i < L; i++) {
    for (j = 1; j < R; j++) {
      if (setup->point_type[i][j] == HVC)
        setup->v[0][i][j] = setup->v[1][i][j] = setup->xtal_HV;
      if (setup->point_type[i][j] == PC)
        setup->v[0][i][j] = setup->v[1][i][j] = 0.0;
    }
  }

  if (!old_setup || !old_setup->fully_depleted) ev_relax_undep(setup);
  else do_relax(setup, 1);
  if (setup->write_field) write_ev(setup);

  if (setup->fully_depleted) {
    printf("Detector is fully depleted.\n");
    /* save potential close to point contact, to use later when calculating depletion voltage */
    for (i = 1; i < L; i++) {
      for (j = 1; j < R; j++) {
        setup->vsave[i][j] = fabs(setup->v[1][i][j]);
      }
    }
  } else {
    printf("Detector is not fully depleted.\n");
    if (setup->bubble_volts > 0)
      printf("Pinch-off bubble at %.1f V potential\n", setup->bubble_volts);
    if (!old_setup) {
      // write a little file that shows any undepleted voxels in the crystal
      FILE *file = fopen("undepleted.txt", "w");
      for (j = R-1; j > 0; j--) {
	setup->undepleted[j][L-1] = '\0';
	fprintf(file, "%s\n", setup->undepleted[j]+1);
      }
      fclose(file);
    }
  }

  return 0;
} /* ev_calc */

/* -------------------------------------- wp_calc ------------------- */
int wp_calc(MJD_Siggen_Setup *setup, MJD_Siggen_Setup *old_setup) {
  int    i, j, pinched_off = 0;
  float  grid = setup->xtal_grid;
  int    L  = lrint(setup->xtal_length/grid)+3;
  int    R  = lrint(setup->xtal_radius/grid)+3;

  if (!old_setup) {
    printf("\n\n ---- starting WP calculation --- \n");
    for (i = 1; i < L; i++) {
      for (j = 1; j < R; j++) {
        setup->v[0][i][j] = setup->v[1][i][j] = 0.01;
      }
    }
  }
  if (old_setup) interpolate(setup, old_setup);

  /* set boundary voltages */
  for (i = 1; i < L; i++) {
    for (j = 1; j < R; j++) {
      if (setup->point_type[i][j] == HVC)
        setup->v[0][i][j] = setup->v[1][i][j] = 0.0;
      else if (setup->point_type[i][j] == PC)
        setup->v[0][i][j] = setup->v[1][i][j] = 1.0;
      else if (setup->undepleted[j][i] == '*') {
        setup->point_type[i][j] = PC;
        setup->v[0][i][j] = setup->v[1][i][j] = 1.0;
      } else if (setup->undepleted[j][i] == 'B') {
        setup->point_type[i][j] = PINCHOFF;
        pinched_off = 1;
      }
    }
  }

  if (pinched_off) wp_relax_undep(setup);
  else do_relax(setup, 0);

  if (setup->write_WP) write_wp(setup);

  return 0;
} /* wp_calc */

/* -------------------------------------- report_config ------------------- */
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
} /* report_config */

/* -------------------------------------- write_ev ------------------- */
int write_ev(MJD_Siggen_Setup *setup) {
  int    i, j, new=1;
  float  grid = setup->xtal_grid;
  int    L  = lrint(setup->xtal_length/grid) + 2;
  int    R  = lrint(setup->xtal_radius/grid) + 2;
  float  r, z;
  float  E_r, E_z, E;
  FILE   *file;
  double ***v = setup->v;

  setup->Emin = 99999.9;
  setup->rmin = setup->zmin = 999.9;

  if (setup->impurity_z0 > 0) {
    // swap voltages back to negative for n-type material
    for (i=1; i<L; i++) {
      for (j=1; j<R; j++) {
        setup->v[new][i][j] = -setup->v[new][i][j];
      }
    }
  }

  /* write potential and field to output file */
  if (!(file = fopen(setup->field_name, "w"))) {
    printf("ERROR: Cannot open file %s for electric field...\n", setup->field_name);
    return 1;
  } else {
    printf("Writing electric field data to file %s\n", setup->field_name);
  }
  /* copy configuration parameters to output file */
  report_config(file, setup->config_file_name);
  fprintf(file, "#\n# HV bias in fieldgen: %.1f V\n", setup->xtal_HV);
  if (setup->fully_depleted) {
    fprintf(file, "# Detector is fully depleted.\n");
  } else {
    fprintf(file, "# Detector is not fully depleted.\n");
    if (setup->bubble_volts > 0.0f)
      fprintf(file, "# Pinch-off bubble at %.0f V potential\n", setup->bubble_volts);
  }
  fprintf(file, "#\n## r (mm), z (mm), V (V),  E (V/cm), E_r (V/cm), E_z (V/cm)\n");
  
  for (j = 1; j < R; j++) {
    r = (j-1) * grid;
    for (i = 1; i < L; i++) {
      z = (i-1) * grid;
      // calc E in r-direction
      if (j == 1) {
        E_r = 0;
      } else if (j == R) {
        E_r = (v[new][i][j-1] - v[new][i][j])/(0.1*grid);
      } else {
        E_r = (v[new][i][j-1] - v[new][i][j+1])/(0.2*grid);
      }
      // calc E in z-direction
      if (i == 1) {
        E_z = (v[new][i][j] - v[new][i+1][j])/(0.1*grid);
      } else if (i == L) {
        E_z = (v[new][i-1][j] - v[new][i][j])/(0.1*grid);
      } else {
        E_z = (v[new][i-1][j] - v[new][i+1][j])/(0.2*grid);
      }
      E = sqrt(E_r*E_r + E_z*E_z);
      fprintf(file, "%7.2f %7.2f %7.1f %7.1f %7.1f %7.1f\n",
              r, z, v[new][i][j], E, E_r, E_z);

      /* check for minimum field inside bulk of detector */
      int k = 3.0/grid;
      if (E > 0.1 && E < setup->Emin &&
          i > k+1 && j < R-k-1 && i < L-k-1 &&
          setup->point_type[i][j] == INSIDE &&
          setup->point_type[i + k][j] == INSIDE &&  // point is at leat 3 mm from a boundary
          setup->point_type[i - k][j] == INSIDE &&
          setup->point_type[i][j + k] == INSIDE &&
          (j < k+1 || setup->point_type[i][j - k] == INSIDE)) {
        setup->Emin = E;
        setup->rmin = r;
        setup->zmin = z;
      }
    }
    fprintf(file, "\n");
  }
  fclose(file);
  printf("\n Minimum bulk field = %.1f V/cm at (r,z) = (%.1f, %.1f) mm\n\n",
         setup->Emin, setup->rmin, setup->zmin);

  if (1) { /* write point_type to output file */
    file = fopen("fields/point_type.dat", "w");
    for (j = 1; j < R; j++) {
      for (i = 1; i < L; i++)
        fprintf(file, "%7.2f %7.2f %2d\n",
                (j-1)*grid, (i-1)*grid, setup->point_type[i][j]);
      fprintf(file, "\n");
    }
    fclose(file);
  }

  return 0;
 } /* write_ev */

/* -------------------------------------- write_wp ------------------- */
 int write_wp(MJD_Siggen_Setup *setup) {
  int    i, j, new=1;;
  float  grid = setup->xtal_grid;
  int    L  = lrint(setup->xtal_length/grid) + 2;
  int    R  = lrint(setup->xtal_radius/grid) + 2;
  float  r, z;
  FILE *file;

  if (!(file = fopen(setup->wp_name, "w"))) {
    printf("ERROR: Cannot open file %s for weighting potential...\n", setup->wp_name);
    return 1;
  } else {
    printf("Writing weighting potential to file %s\n\n", setup->wp_name);
  }

  /* copy configuration parameters to output file */
  report_config(file, setup->config_file_name);
  fprintf(file, "#\n# HV bias in fieldgen: %.1f V\n", setup->xtal_HV);
  if (setup->fully_depleted) {
    fprintf(file, "# Detector is fully depleted.\n");
  } else {
    fprintf(file, "# Detector is not fully depleted.\n");
    if (setup->bubble_volts > 0.0f)
      fprintf(file, "# Pinch-off bubble at %.0f V potential\n", setup->bubble_volts);
  }
  fprintf(file, "#\n## r (mm), z (mm), WP\n");
  for (j = 1; j < R; j++) {
    r = (j-1) * grid;
    for (i = 1; i < L; i++) {
      z = (i-1) * grid;
      fprintf(file, "%7.2f %7.2f %12.6e\n", r, z, setup->v[new][i][j]);
    }
    fprintf(file, "\n");
  }
  fclose(file);

  return 0;
 } /* write_wp */

/* -------------------------------------- grid_init ------------------- */
#define SQ(x) ((x)*(x))
int grid_init(MJD_Siggen_Setup *setup) {
  float  grid = setup->xtal_grid;
  int    L  = lrint(setup->xtal_length/grid)+3;
  int    R  = lrint(setup->xtal_radius/grid)+3;
  int    i, j;
  float  r, z;

  /* first malloc arrays in setup */
  if ((setup->impurity   = malloc(L * sizeof(*setup->impurity)))   == NULL ||
      (setup->eps        = malloc(L * sizeof(*setup->eps)))        == NULL ||
      (setup->eps_dr     = malloc(L * sizeof(*setup->eps_dr)))     == NULL ||
      (setup->eps_dz     = malloc(L * sizeof(*setup->eps_dz)))     == NULL ||
      (setup->v[0]       = malloc(L * sizeof(*setup->v[0])))       == NULL ||
      (setup->v[1]       = malloc(L * sizeof(*setup->v[1])))       == NULL ||
      (setup->undepleted = malloc(R * sizeof(*setup->undepleted))) == NULL ||
      (setup->s1         = malloc(R * sizeof(*setup->s1)))         == NULL ||
      (setup->s2         = malloc(R * sizeof(*setup->s2)))         == NULL ||
      (setup->vsave      = malloc(L * sizeof(*setup->vsave)))      == NULL ||
      (setup->point_type = malloc(L * sizeof(*setup->point_type))) == NULL) {
    printf("malloc failed\n");
    return -1;
  }
  /* start from i=1 so that i=0 can be used for reflection symmetry around r=0 or z=0 */
  for (i = 1; i < L; i++) {
    if ((setup->impurity[i]   = malloc(R * sizeof(**setup->impurity)))   == NULL ||
        (setup->v[0][i]       = malloc(R * sizeof(**setup->v[0])))       == NULL ||
	(setup->v[1][i]       = malloc(R * sizeof(**setup->v[1])))       == NULL ||
	(setup->vsave[i]      = malloc(R * sizeof(**setup->vsave)))      == NULL ||
	(setup->point_type[i] = malloc(R * sizeof(**setup->point_type))) == NULL) {
      printf("malloc failed\n");
      return -1;
    }
  }
  for (i = 1; i < L; i++) {
    if ((setup->eps[i]        = malloc(R * sizeof(**setup->eps)))    == NULL ||
        (setup->eps_dr[i]     = malloc(R * sizeof(**setup->eps_dr))) == NULL ||
        (setup->eps_dz[i]     = malloc(R * sizeof(**setup->eps_dz))) == NULL) {
      printf("malloc failed\n");
      return -1;
    }
    for (j = 0; j < R; j++)
      setup->eps[i][j] = setup->eps_dz[i][j] = setup->eps_dr[i][j] = 16.0;
  }
  for (j = 0; j < R; j++) {
    if ((setup->undepleted[j] = malloc(L * sizeof(**setup->undepleted))) == NULL) {
      printf("malloc failed\n");
      return -1;
    }
    memset(setup->undepleted[j], ' ', L);
  }
        
  /* set up reflection symmetry around r=0 or z=0 */
  setup->impurity[0]   = malloc(R * sizeof(**setup->impurity));
  setup->v[0][0]   = setup->v[0][2];
  setup->v[1][0]   = setup->v[1][2];
  setup->eps[0]    = setup->eps[1];
  setup->eps_dr[0] = setup->eps_dr[1];
  setup->eps_dz[0] = setup->eps_dz[1];
  setup->point_type[0] = setup->point_type[1];

  /* ------------------------------------------------------------ */
  /* weighting values for the relaxation alg. as a function of r
     in the following we divide areas and volumes by pi
     r_bin   rmax  A_top A_outside A_inside  volume  total_surf  out/top  tot/vol
     0     1/2    1/4      1         0       1/4      1.5         4        6  << special case
     1     3/2      2      3         1        2        8        3/2        4
     2     5/2      4      5         3        4       16        5/4        4
     3     7/2      6      7         5        6       24        7/6        4
     r   r+0.5     2r    2r+1      2r-1      2r       8r     (2r+1)/2r     4
     = 1+0.5/r
  */
  setup->s1[1] = 4.0;
  setup->s2[1] = 0.0;
  for (i=2; i<R; i++) {
    setup->s1[i] = 1.0 + 0.5 / (double) (i-1);   //  for r+1
    setup->s2[i] = 1.0 - 0.5 / (double) (i-1);   //  for r-1
  }
  setup->s2[1] = setup->s1[1]; // special case for reflection symm at r=0

  /* set up point types for boundary conditions etc */
  float lith  = setup->Li_thickness;
  float zmax  = setup->xtal_length - lith;
  float rmax  = setup->xtal_radius - lith;
  float r1, z1, br, a, b;

  for (i = 1; i < L; i++) {
    z = (i-1) * grid;
    z1 = zmax - z;  // distance from top of crystal
    for (j = 1; j < R; j++) {
      r = (j-1) * grid;
      r1 = rmax - r;  // distance from outer radius

      setup->point_type[i][j] = INSIDE;

      // outside (HV) contact:
      br = setup->hole_bullet_radius + 0.8*lith;
      a = setup->hole_radius + lith - br;
      b = zmax - setup->hole_length + br;
      if (z >= zmax || r >= rmax ||
          (r >= setup->wrap_around_radius && z <= lith) ||      // wrap-around  // COULD BE r >
          /* check 45-degree bottom outer taper of crystal */
          (z  <= setup->bottom_taper_length + 0.71*lith &&
           r1 <= setup->bottom_taper_length + 0.71*lith - z) ||
          /* check hole */
          (r  <= setup->hole_radius + lith &&
           z1 <= setup->hole_length &&  // note no lith added here, since it was subtracted from zmax
            /* check hole bulletization */
            (r <= a || z >= b || SQ(b-z) + SQ(a-r) <= SQ(br))) ||
          /* check inner taper of hole */
          (z1 <= setup->inner_taper_length &&
           r  <= setup->hole_radius + lith +
                  ((setup->inner_taper_length - z1) *
                   setup->inner_taper_width / setup->inner_taper_length)) ||
          /* check outer top taper of crystal */
          (z1 <= setup->outer_taper_length &&
           r1 <= ((setup->outer_taper_length - z1) *
                   setup->outer_taper_width / setup->outer_taper_length)) ||
          /* check 45-degree bottom outer taper of crystal */
          (z  <= setup->bottom_taper_length + 0.71*lith &&
           r1 <= setup->bottom_taper_length + 0.71*lith - z))
        setup->point_type[i][j] = HVC;

      /* check top and bottom bulletizations */
      br = setup->top_bullet_radius - lith * 0.7;
      if (z1 <= br &&r1 <= br &&
          SQ(br - r1) + SQ(br - z1) >= br*br)
        setup->point_type[i][j] = HVC;
      br = setup->bottom_bullet_radius - lith * 0.7;
      if (z <= br + lith && r1 <= br &&
          SQ(br - r1) + SQ(br - z + lith) >= br*br)
        setup->point_type[i][j] = HVC;

      if (setup->point_type[i][j] != INSIDE) continue;

      /* check for inside (point) contact, with optional bulletization */
      if (z <= setup->pc_length && r <= setup->pc_radius) {
        if (setup->bulletize_PC) {
          br = setup->pc_radius;
          // if length <= radius, use length as bulletization radius
          if (setup->pc_length <= setup->pc_radius) br = setup->pc_length;
          a = setup->pc_radius - br;
          b = setup->pc_length - br;
          if (r < a || z < b || SQ(b-z) + SQ(a-r) <= SQ(br)) setup->point_type[i][j] = PC;
        }
        else {
          setup->point_type[i][j] = PC;
        }
      }
      /* ---------------------------------------------- 
         edges of inside contact; if the PC radius and/or length is not in the middle
           of a pixel, we want to modify interpolation of V in surrounding pixels...
           in this case, the radius/length > grid point, so it modifies the
           interpolation for the next point out
           // FIXME: Check for adjacent ditch *//*
        if (r == rrc[z] && drrc[z] < -0.05) {
          bulk[z][r] = 1;  // flag for radial edge of PC
          frrc[z] = -1.0/drrc[z];  // interpolation weight for pixel at (r-1)
          // only part of the pixel has volume charge density, the rest is contact
          vfraction[z][r] *= -2.0*drrc[z];
        }
        // z edge of inside contact; if the PC length is not in the middle
        //   of a pixel, we want to modify interpolation of V in surrounding pixels
        if (z == LC && dLC < -0.05) {
          bulk[z][r] = 2;  // flag for z edge of PC
          fLC = -1.0/dLC;  // interpolation weight for pixel at (z-1)
          // only part of the pixel has volume charge density, the rest is contact
          vfraction[z][r] *= -2.0*dLC;
        }
        ---------------------------------------------- */

      /* check for inside ditch
         boundary condition at Ge-vacuum interface:
         epsilon0 * E_vac = espilon_Ge * E_Ge
      */
      if (setup->ditch_depth  > 0 && z <= setup->ditch_depth  &&     // COULD BE z <
          r < setup->wrap_around_radius &&                           // COULD BE r <=
          r > setup->wrap_around_radius - setup->ditch_thickness) {  // COULD BE r >=
        setup->point_type[i][j] = DITCH;
        setup->eps[i][j] = setup->eps_dz[i][j] = setup->eps_dr[i][j] = 1.0;
      }

    }
  }

  /* for pixels adjacent to the ditch, set point_type to DITCH_EDGE
     andfor z=0, set flag for passivated surface */
  for (i = 1; i < L; i++) {
    for (j = 1; j < R; j++) {
      setup->eps_dr[i][j-1] = (setup->eps[i][j-1] + setup->eps[i][j]) / 2.0f;
      setup->eps_dz[i-1][j] = (setup->eps[i-1][j] + setup->eps[i][j]) / 2.0f;
      if (setup->point_type[i][j] == INSIDE &&
          (setup->point_type[i-1][j] == DITCH ||
           setup->point_type[i][j-1] == DITCH ||
           setup->point_type[i][j+1] == DITCH)) setup->point_type[i][j] = DITCH_EDGE;
      if (i == 1 && setup->point_type[i][j] == INSIDE &&
          (j-1) * grid < setup->wrap_around_radius) setup->point_type[i][j] = PASSIVE;
    }
    setup->eps_dr[i][0] = setup->eps_dr[i][1];
  }

  /* set up impurity array */
  double *imp_z, imp_ra = 0, imp_rm = 1;
  double e_over_E = 11.310; // e/epsilon; for 1 mm2, charge units 1e10 e/cm3, espilon = 16*epsilon0
  /* malloc local array */
  if ((imp_z  = malloc(L * sizeof(*imp_z))) == NULL) {
    printf("malloc failed\n");
    return -1;
  }

  if (setup->rho_z_spe[0] == 0) {
    for (i = 1; i < L; i++) {
      z = (i-1) * grid;
      imp_z[i] = e_over_E * grid*grid / 4.0 *
                 (setup->impurity_z0 +
                  setup->impurity_gradient * z * 0.1 +
                  setup->impurity_quadratic *
                  (1.0 - SQ(z - setup->xtal_length/2.0) / SQ(setup->xtal_length/2.0)));
    }
  } else {
    for (i = 1; i < L; i++)
      imp_z[i] = e_over_E * grid*grid / 4.0 * setup->rho_z_spe[(int) ((i-1) * grid)];
  }
  if (0) printf("imp_z: %.2f %.2f %.2f\n",
                 imp_z[1]/(e_over_E * grid*grid / 4.0),
                 imp_z[L/2]/(e_over_E * grid*grid / 4.0),
                 imp_z[L-2]/(e_over_E * grid*grid / 4.0));
  for (j = 1; j < R; j++) {
    r = (j-1) * grid;
    if (setup->impurity_rpower > 0.1) {
      imp_ra = setup->impurity_radial_add * e_over_E *
        pow((double) r / setup->xtal_radius, setup->impurity_rpower);
      imp_rm = 1.0 + (setup->impurity_radial_mult - 1.0f) *
        pow((double) r / setup->xtal_radius, setup->impurity_rpower);
    }
    for (i = 1; i < L; i++)  setup->impurity[i][j] = imp_z[i] * imp_rm + imp_ra;
    if (setup->point_type[1][j] == PASSIVE) {
      setup->impurity[1][j] += setup->impurity_surface * e_over_E * grid;
      // special hack for Ponama 1 in GALATEA?
      // if (r < 8) setup->impurity[1][j] -= 2.0 * e_over_E * grid;
      // setup->impurity[1][j] -= 0.3 * 13.0/r * e_over_E * grid;
      // setup->impurity[1][j] -= 0.3 * (1.0 - r/15.0) * e_over_E * grid;
    }
  }

  /* free local and no-longer-needed arrays */
  free(imp_z);
  for (i = 1; i < L; i++) free(setup->eps[i]);
  free(setup->eps);

  return 0;   
} /* grid_init */
#undef SQ

/* -------------------------------------- do_relax ------------------- */
int do_relax(MJD_Siggen_Setup *setup, int ev_calc) {
  int    old = 1, new = 0, iter, r, z;
  float  grid = setup->xtal_grid;
  int    L  = lrint(setup->xtal_length/grid)+1;
  int    R  = lrint(setup->xtal_radius/grid)+1;
  double eps_sum, v_sum, mean, save_dif;
  double dif, sum_dif, max_dif;
  double ***v = setup->v, **eps_dr = setup->eps_dr, **eps_dz = setup->eps_dz;
  double *s1 = setup->s1, *s2 = setup->s2;
  double e_over_E = 11.310; // e/epsilon; for 1 mm2, charge units 1e10 e/cm3, espilon = 16*epsilon0


  if (setup->vacuum_gap > 0) {   // save impurity value along passivated surface
    for (r = 1; r < R; r++)
      setup->impurity[0][r] = setup->impurity[1][r];
  }

  for (iter = 0; iter < setup->max_iterations; iter++) {

    // the following definition of the factor for over-relaxation improves convergence
    //     time by a factor ~ 70-120 for a 2kg ICPC detector, grid = 0.1 mm
    //   OR_fact increases with increasing volxel count (L*R)
    //         and with increasing iteration number
    //   0.997 is maximum asymptote for very large pixel count and iteration number
    double OR_fact = ((0.997 - 600.0/(L*R)) * (1.0 - 0.9/(double)(1+iter/6)));
    if (600.0/(L*R) > 0.5)
      OR_fact = (0.5 * (1.0 - 0.9/(double)(1+iter/6)));
    if (iter < 2) OR_fact = 0.0;

    old = new;
    new = 1 - new;
    sum_dif = 0;
    max_dif = 0;

    if (setup->vacuum_gap > 0) {   // modify impurity value along passivated surface
      for (r = 1; r < R; r++)      //   due to surface charge induced by capacitance
        if (ev_calc)
          setup->impurity[1][r] = setup->impurity[0][r] +
            v[old][1][r] * 5.52e-4 * e_over_E * grid / setup->vacuum_gap;
        else 
          setup->impurity[1][r] = v[old][1][r] * 5.52e-4 * e_over_E * grid / setup->vacuum_gap;
    }

    /* start from z=1 and r=1 so that (z,r)=0 can be
       used for reflection symmetry around r=0 or z=0 */
    for (z = 1; z < L; z++) {
      /* manage r=0 reflection symmetry */
      setup->v[old][z][0] = setup->v[old][z][2];

      for (r = 1; r < R; r++) {
        if (setup->point_type[z][r] < INSIDE) continue;   // HV or point contact
        save_dif = v[old][z][r] - v[new][z][r];      // step difference from previous iteration

        if (setup->point_type[z][r] < DITCH) {       // normal bulk or passivated surface, no complications
          v_sum = (v[old][z+1][r] + v[old][z][r+1]*s1[r] +
                   v[old][z-1][r] + v[old][z][r-1]*s2[r]);
          eps_sum = 4;
          if (r == 1) eps_sum = 2 + s1[r] + s2[r];
        } else if (setup->point_type[z][r] >= DITCH) {  // in or adjacent to the ditch
          v_sum = (v[old][z+1][r]*eps_dz[z  ][r] + v[old][z][r+1]*eps_dr[z][r  ]*s1[r] +
                   v[old][z-1][r]*eps_dz[z-1][r] + v[old][z][r-1]*eps_dr[z][r-1]*s2[r]);
          eps_sum = (eps_dz[z][r]   + eps_dr[z][r]  *s1[r] +
                     eps_dz[z-1][r] + eps_dr[z][r-1]*s2[r]);
          /*
        } else if (setup->point_type[z][r] == 1) {    // interpolated radial edge of point contact
          // since the PC radius is not in the middle of a pixel,
          //   use a modified weight for the interpolation to (r-1)
          v_sum = v[old][z+1][r]*eps_dz[z][r] + v[old][z][r+1]*eps_dr[z][r]*s1[r] +
            v[old][z][r-1]*eps_dr[z][r-1]*s2[r];// *frrc[z];
          eps_sum = eps_dz[z][r] + eps_dr[z][r]*s1[r] + eps_dr[z][r-1]*s2[r];// *frrc[z];
          if (z > 0) {
            v_sum += v[old][z-1][r]*eps_dz[z-1][r];
            eps_sum += eps_dz[z-1][r];
          } else {
            v_sum += v[old][z+1][r]*eps_dz[z][r];  // reflection symm around z=0
            eps_sum += eps_dz[z][r];
          }
        } else if (setup->point_type[z][r] == 2) {    // interpolated z edge of point contact
          // since the PC length is not in the middle of a pixel,
          //   use a modified weight for the interpolation to (z-1)
          v_sum = v[old][z+1][r]*eps_dz[z][r] + v[old][z][r+1]*eps_dr[z][r]*s1[r] +
            v[old][z-1][r]*eps_dz[z-1][r];// *fLC;
          eps_sum = eps_dz[z][r] + eps_dr[z][r]*s1[r] + eps_dz[z-1][r];// *fLC;
          if (r > 0) {
            v_sum += v[old][z][r-1]*eps_dr[z][r-1]*s2[r];
            eps_sum += eps_dr[z][r-1]*s2[r];
          } else {
            v_sum += v[old][z][r+1]*eps_dr[z][r]*s1[r];  // reflection symm around r=0
            eps_sum += eps_dr[z][r]*s1[r];
          }
          // check for cases where the PC corner needs modification in both r and z
          if (z == LC && setup->point_type[z-1][r] == 1) {
            v_sum += v[old][z][r-1]*eps_dr[z][r-1]*s2[r];// *(frrc[z]-1.0);
            eps_sum += eps_dr[z][r-1]*s2[r];// *(frrc[z]-1.0);
          }
          */
        }

        // calculate the interpolated mean potential and the effect of the space charge
        mean = v_sum / eps_sum;
        if (ev_calc ||
            (setup->vacuum_gap > 0 && z == 1))
          v[new][z][r] = mean + setup->impurity[z][r];
        else
          v[new][z][r] = mean;

        // calculate difference from last iteration, for convergence check
        dif = v[old][z][r] - v[new][z][r];
        if (dif < 0) dif = -dif;
        sum_dif += dif;
        if (max_dif < dif) max_dif = dif;
        // do over-relaxation
        v[new][z][r] += OR_fact*save_dif;

        //if (v[new][z][r] < 0) v[new][z][r] = 0;
      }
    }

    // report results for some iterations
    if (iter < 10 || (iter < 600 && iter%100 == 0) || iter%1000 == 0) {
      if (0 && ev_calc) {
        printf("%5d %d %d %.10f %.10f\n", iter, old, new, max_dif, sum_dif/(L-2)/(R-2));
      } else {
        printf("%5d %d %d %.10f %.10f ; %.10f %.10f\n",
               iter, old, new, max_dif, sum_dif/(L-2)/(R-2),
               v[new][L/2][R/2], v[new][L/3][R/3]);
      }
    }
    // check for convergence
    if ( ev_calc && max_dif < 0.00000008) break;
    if (!ev_calc && max_dif < 0.0000000001) break;

  }

  printf(">> %d %.16f\n\n", iter, sum_dif);
  if (setup->vacuum_gap > 0) {   // restore impurity value along passivated surface
    for (r = 1; r < R; r++)
      setup->impurity[1][r] = setup->impurity[0][r];
  }

  return 0;
} /* do_relax */

/* -------------------------------------- ev_relax_undep ------------------- */
int ev_relax_undep(MJD_Siggen_Setup *setup) {
  int    old = 1, new = 0, iter, r, z;
  float  grid = setup->xtal_grid;
  int    L  = lrint(setup->xtal_length/grid)+1;
  int    R  = lrint(setup->xtal_radius/grid)+1;
  double eps_sum, v_sum, mean, save_dif, min;
  double dif, sum_dif, max_dif, bubble_volts, bv2, bv3;
  double ***v = setup->v, **eps_dr = setup->eps_dr, **eps_dz = setup->eps_dz;
  double *s1 = setup->s1, *s2 = setup->s2;
  char   **undep = setup->undepleted;
  double e_over_E = 11.310; // e/epsilon; for 1 mm2, charge units 1e10 e/cm3, espilon = 16*epsilon0


  if (setup->vacuum_gap > 0) {   // save impurity value along passivated surface
    for (r = 1; r < R; r++)
      setup->impurity[0][r] = setup->impurity[1][r];
  }

#ifndef BUBBLE1
  /* initialise the undepleted array for use with bubble depletion */
  for (z = 1; z < L; z++) {
    for (r = 1; r < R; r++) {
      if (setup->point_type[z][r] >= INSIDE) undep[r][z] = 0;
    }
  }
#endif

  for (iter = 0; iter < setup->max_iterations; iter++) {

    double OR_fact = ((0.997 - 600.0/(L*R)) * (1.0 - 0.9/(double)(1+iter/6)));
    if (600.0/(L*R) > 0.5)
      OR_fact = (0.5 * (1.0 - 0.9/(double)(1+iter/6)));
    if (iter < 2) OR_fact = 0.0;

    old = new;
    new = 1 - new;
    sum_dif = 0;
    max_dif = 0;
    bubble_volts = 0;
    //if (bv3 > 0) bubble_volts = bv2/bv3 + 0.01;
    bv2 = 0;
    bv3 = 0;

    if (setup->vacuum_gap > 0) {   // modify impurity value along passivated surface
      for (r = 1; r < R; r++)      //   due to surface charge induced by capacitance
        setup->impurity[1][r] = setup->impurity[0][r] +
          v[old][1][r] * 5.52e-4 * e_over_E * grid / setup->vacuum_gap;
    }

    /* start from z=1 and r=1 so that (z,r)=0 can be
       used for reflection symmetry around r=0 or z=0 */
    for (z = 1; z < L; z++) {
      /* manage r=0 reflection symmetry */
      setup->v[old][z][0] = setup->v[old][z][2];

      for (r = 1; r < R; r++) {
        if (setup->point_type[z][r] < INSIDE) continue;   // HV or point contact
        save_dif = v[old][z][r] - v[new][z][r];      // step difference from previous iteration

        if (setup->point_type[z][r] < DITCH) {       // normal bulk or passivated surface, no complications
          v_sum = (v[old][z+1][r] + v[old][z][r+1]*s1[r] +
                   v[old][z-1][r] + v[old][z][r-1]*s2[r]);
          eps_sum = 4;
          if (r == 1) eps_sum = 2 + s1[r] + s2[r];
        } else if (setup->point_type[z][r] >= DITCH) {  // in or adjacent to the ditch
          v_sum = (v[old][z+1][r]*eps_dz[z  ][r] + v[old][z][r+1]*eps_dr[z][r  ]*s1[r] +
                   v[old][z-1][r]*eps_dz[z-1][r] + v[old][z][r-1]*eps_dr[z][r-1]*s2[r]);
          eps_sum = (eps_dz[z][r]   + eps_dr[z][r]  *s1[r] +
                     eps_dz[z-1][r] + eps_dr[z][r-1]*s2[r]);
        }

        // calculate the interpolated mean potential and the effect of the space charge
        min = fminf(fminf(v[old][z+1][r], v[old][z][r+1]),
                    fminf(v[old][z-1][r], v[old][z][r-1]));
        mean = v_sum / eps_sum;
        v[new][z][r] = mean + setup->impurity[z][r];

#ifndef BUBBLE1
        undep[r][z] /= 2;
        if (v[new][z][r] <= 0) {
          v[new][z][r] = 0;
          undep[r][z] = 4;  // do not do over-relaxation for 3 iterations
        } else if (v[new][z][r] <= min) {
          if (bubble_volts == 0) bubble_volts = min + 0.01;
          v[new][z][r] = bubble_volts;
          bv2 += min;
          bv3 += 1.0;
          undep[r][z] = 8;  // do not do over-relaxation for 4 iterations
        }
#else
	undep[r][z] = '.';
        if (v[new][z][r] <= 0) {
          v[new][z][r] = 0;
          undep[r][z] = '*';
          save_dif = 0;  // do not do over-relaxation
        } else if (v[new][z][r] < min) {
          if (bubble_volts == 0) bubble_volts = min + 0.01;
          v[new][z][r] = bubble_volts;
          undep[r][z] = '*';
          save_dif = 0;  // do not do over-relaxation
        }
#endif

        // calculate difference from last iteration, for convergence check
        dif = v[old][z][r] - v[new][z][r];
        if (dif < 0) dif = -dif;
        sum_dif += dif;
        if (max_dif < dif) max_dif = dif;
        // do over-relaxation
#ifndef BUBBLE1
        if (!undep[r][z])
#endif
        v[new][z][r] += OR_fact*save_dif;
      }
    }

    // report results for some iterations
    if (iter < 10 || (iter < 600 && iter%100 == 0) || iter%1000 == 0) {
      if (0) {
        printf("%5d %d %d %.10f %.10f\n", iter, old, new, max_dif, sum_dif/(L-2)/(R-2));
      } else {
        printf("%5d %d %d %.10f %.10f ; %.10f %.10f bubble %.2f %.0f\n",
               iter, old, new, max_dif, sum_dif/(L-2)/(R-2),
               v[new][L/2][R/2], v[new][L/3][R/3], bubble_volts, bv3);
      }
    }
    // check for convergence
    if (max_dif < 0.00000008) break;
 
  }
  printf(">> %d %.16f\n\n", iter, sum_dif);

  setup->bubble_volts = bubble_volts;
  setup->fully_depleted = 1;
  for (r=1; r<R; r++) {
    for (z=1; z<L; z++) {
#ifndef BUBBLE1
      if (setup->point_type[z][r] < INSIDE) {
        undep[r][z] = ' ';
      } else if (undep[r][z] == 0) {
        undep[r][z] = '.';
      } else {
        if (undep[r][z] > 4) undep[r][z] = 'B';  // identifies pinch-off
        else undep[r][z] = '*';
        setup->fully_depleted = 0;
      }
#else
     if (undep[r][z] == '*') {
        setup->fully_depleted = 0;
        if (v[new][z][r] > 0.001) undep[r][z] = 'B';  // identifies pinch-off
      }
#endif
    }
  }

  printf("bubble %.1f\n", bubble_volts);
  if (setup->vacuum_gap > 0) {   // restore impurity value along passivated surface
    for (r = 1; r < R; r++)
      setup->impurity[1][r] = setup->impurity[0][r];
  }

  return 0;
} /* ev_relax_undep */

/* -------------------------------------- wp_relax_undep ------------------- */
int wp_relax_undep(MJD_Siggen_Setup *setup) {
  int    old = 1, new = 0, iter, r, z;
  float  grid = setup->xtal_grid;
  int    L  = lrint(setup->xtal_length/grid)+1;
  int    R  = lrint(setup->xtal_radius/grid)+1;
  double eps_sum, v_sum, mean, save_dif, pinched_sum1, pinched_sum2;
  double dif, sum_dif, max_dif;
  double ***v = setup->v, **eps_dr = setup->eps_dr, **eps_dz = setup->eps_dz;
  double *s1 = setup->s1, *s2 = setup->s2;
  double e_over_E = 11.310; // e/epsilon; for 1 mm2, charge units 1e10 e/cm3, espilon = 16*epsilon0


  if (setup->vacuum_gap > 0) {   // save impurity value along passivated surface
    for (r = 1; r < R; r++)
      setup->impurity[0][r] = setup->impurity[1][r];
  }

  for (iter = 0; iter < setup->max_iterations; iter++) {

   double OR_fact = ((0.997 - 600.0/(L*R)) * (1.0 - 0.9/(double)(1+iter/6)));
    if (600.0/(L*R) > 0.5)
      OR_fact = (0.5 * (1.0 - 0.9/(double)(1+iter/6)));
    if (iter < 2) OR_fact = 0.0;

    old = new;
    new = 1 - new;
    sum_dif = 0;
    max_dif = 0;
    pinched_sum1 = pinched_sum2 = 0.0;

    if (setup->vacuum_gap > 0) {   // modify impurity value along passivated surface
      for (r = 1; r < R; r++)      //   due to surface charge induced by capacitance
        setup->impurity[1][r] = v[old][1][r] * 5.52e-4 * e_over_E * grid / setup->vacuum_gap;
    }

    /* start from z=1 and r=1 so that (z,r)=0 can be
       used for reflection symmetry around r=0 or z=0 */
    for (z = 1; z < L; z++) {
      /* manage r=0 reflection symmetry */
      setup->v[old][z][0] = setup->v[old][z][2];

      for (r = 1; r < R; r++) {
        if (setup->point_type[z][r] < INSIDE) continue;   // HV or point contact
        save_dif = v[old][z][r] - v[new][z][r];      // step difference from previous iteration

        if (setup->point_type[z][r] < PINCHOFF) {       // normal bulk or passivated surface, no complications
          v_sum = (v[old][z+1][r] + v[old][z][r+1]*s1[r] +
                   v[old][z-1][r] + v[old][z][r-1]*s2[r]);
          eps_sum = 4;
          if (r == 1) eps_sum = 2 + s1[r] + s2[r];
        } else if (setup->point_type[z][r] == PINCHOFF) {  // in or adjacent to the ditch
          if (setup->point_type[z+1][r] < PINCHOFF) {
            pinched_sum1 += v[old][z+1][r]*eps_dz[z][r];
            pinched_sum2 += eps_dz[z][r];
          }
          if (setup->point_type[z][r+1] < PINCHOFF) {
            pinched_sum1 += v[old][z][r+1]*eps_dr[z][r]*s1[r];
            pinched_sum2 += eps_dr[z][r]*s1[r];
          }
          if (setup->point_type[z-1][r] < PINCHOFF) {
            pinched_sum1 += v[old][z-1][r]*eps_dz[z-1][r];
            pinched_sum2 += eps_dz[z-1][r];
          }
          if (setup->point_type[z][r-1] < PINCHOFF) {
            pinched_sum1 += v[old][z][r-1]*eps_dr[z][r-1]*s2[r];
            pinched_sum2 += eps_dr[z][r-1]*s2[r];
          }
          v_sum = pinched_sum1;
          eps_sum = pinched_sum2;
        } else if (setup->point_type[z][r] >= DITCH) {  // in or adjacent to the ditch
          v_sum = (v[old][z+1][r]*eps_dz[z  ][r] + v[old][z][r+1]*eps_dr[z][r  ]*s1[r] +
                   v[old][z-1][r]*eps_dz[z-1][r] + v[old][z][r-1]*eps_dr[z][r-1]*s2[r]);
          eps_sum = (eps_dz[z][r]   + eps_dr[z][r]  *s1[r] +
                     eps_dz[z-1][r] + eps_dr[z][r-1]*s2[r]);
        }

        if (setup->point_type[z][r] != PINCHOFF) {
          // calculate the interpolated mean potential and the effect of the space charge
          mean = v_sum / eps_sum;
          if (setup->vacuum_gap > 0 && z == 1)
            v[new][z][r] = mean + setup->impurity[z][r];
          else
            v[new][z][r] = mean;

          // calculate difference from last iteration, for convergence check
          dif = v[old][z][r] - v[new][z][r];
          if (dif < 0) dif = -dif;
          sum_dif += dif;
          if (max_dif < dif) max_dif = dif;
          // do over-relaxation
          v[new][z][r] += OR_fact*save_dif;
        }
      }
    }
    if (pinched_sum2 > 0.1) {
      mean = pinched_sum1 / pinched_sum2;
      for (z=1; z<L; z++) {
        for (r=1; r<R; r++) {
          if (setup->point_type[z][r] == PINCHOFF) {
            v[new][z][r] = mean;
            dif = v[old][z][r] - v[new][z][r];
            if (dif < 0) dif = -dif;
            sum_dif += dif;
            if (max_dif < dif) max_dif = dif;
          }
        }
      }
    }

    // report results for some iterations
    if (iter < 10 || (iter < 600 && iter%100 == 0) || iter%1000 == 0) {
      printf("%5d %d %d %.10f %.10f ; %.10f %.10f\n",
             iter, old, new, max_dif, sum_dif/(L-2)/(R-2),
             v[new][L/2][R/2], v[new][L/3][R/3]);
    }
    // check for convergence
    if (max_dif < 0.0000000001) break;

  }

  printf(">> %d %.16f\n\n", iter, sum_dif);
  if (setup->vacuum_gap > 0) {   // restore impurity value along passivated surface
    for (r = 1; r < R; r++)
      setup->impurity[1][r] = setup->impurity[0][r];
  }

  return 0;
} /* wp_relax_undep */

/* -------------------------------------- interpolate ------------------- */
int interpolate(MJD_Siggen_Setup *setup, MJD_Siggen_Setup *old_setup) {
  int    n, i, j, i2, j2, zmin, rmin, zmax, rmax;
  int    L  = lrint(old_setup->xtal_length/old_setup->xtal_grid)+3;
  int    R  = lrint(old_setup->xtal_radius/old_setup->xtal_grid)+3;
  int    L2 = lrint(setup->xtal_length/setup->xtal_grid)+3;
  int    R2 = lrint(setup->xtal_radius/setup->xtal_grid)+3;
  float  f, f1r, f1z, f2r, f2z;
  double ***v = setup->v, **ov = old_setup->v[1];

  /* the previous calculation was on a coarser grid...
     now copy/expand the potential to the new finer grid
  */
  n = (int) (old_setup->xtal_grid / setup->xtal_grid + 0.5);
  f = 1.0 / (float) n;
  printf("\ngrid %.4f -> %.4f; ratio = %d %.3f\n\n",
         old_setup->xtal_grid, setup->xtal_grid, n, f);
  for (i = i2 = 1; i < L-1; i++) {
    zmin = i2;
    zmax = i2 + n;
    if (zmax > L2-1) zmax = L2-1;
    for (j = j2 = 1; j < R-1; j++) {
      f1z = 0.0;
      rmin = j2;
      rmax = j2 + n;
      if (rmax > R2-1) rmax = R2-1;
      for (i2 = zmin; i2 < zmax; i2++) {
        f2z = 1.0 - f1z;
        f1r = 0.0;
        for (j2 = rmin; j2 < rmax; j2++) {
          f2r = 1.0 - f1r;
          v[0][i2][j2] = v[1][i2][j2] =      // linear interpolation
            f2z*f2r*ov[i][j  ] + f1z*f2r*ov[i+1][j  ] +
            f2z*f1r*ov[i][j+1] + f1z*f1r*ov[i+1][j+1];
          f1r += f;
        }
        f1z += f;
      }
      j2 = rmax;
    }
    i2 = zmax;
  }

  return 0;
} /* interpolate */
