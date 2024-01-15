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
      - added -r option to allow reading of z-impurity profile from a .spe file

   Nov 2019
      - major rewrite and re-arragement of the code
      - broke out functions grid_init, ev_calc, wp_calc, write_ev, write_wp, do_relax, interpolate
      - got an additional factor of ~3 speedup
      - more modular and readable
   Dec 2019
      - added code to properly handle pixels next to p+ and n+ contacts
        - much better calculation of V and E near contacts

      - todo now:  execute WD=1 option
*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "mjd_siggen.h"
#include "detector_geometry.h"

#define MAX_ITS 50000     // default max number of iterations for relaxation

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
  int   WD = 0;  // 0: do not write out depletion surface
                 // 1: write out depletion surface to depl_<HV>.dat

  int   i, j, k;
  FILE  *fp;


  if (argc < 2 || argc%2 != 0 || read_config(argv[1], &setup)) {
    printf("Usage: %s <config_file_name> [options]\n"
           "   Possible options:\n"
	   "      -b bias_volts\n"
	   "      -w {0,1}  (do_not/do write the field file)\n"
	   "      -d {0,1}  (do_not/do write the depletion surface)\n"
	   "      -p {0,1}  (do_not/do write the WP file)\n"
           "      -r rho_spectrum_file_name\n"
           "      -z <rho_spectrum_offset_mm>\n", argv[0]);
    return 1;
  }
  strncpy(setup.config_file_name, argv[1], sizeof(setup.config_file_name));

  if (setup.xtal_grid < 0.001) setup.xtal_grid = 0.5;
  BV = setup.xtal_HV;
  WV = setup.write_field;
  setup.rho_z_spe[0] = 0;

  for (i=2; i<argc-1; i++) {
    if (strstr(argv[i], "-b")) {
      BV = setup.xtal_HV = atof(argv[++i]);   // bias volts
    } else if (strstr(argv[i], "-w")) {
      WV = atoi(argv[++i]);               // write-out options
    } else if (strstr(argv[i], "-d")) {
      WD = atoi(argv[++i]);               // write-out options
    } else if (strstr(argv[i], "-p")) {
      setup.write_WP = atoi(argv[++i]);   // weighting-potential options
    } else if (strstr(argv[i], "-r")) {
      if (!(fp = fopen(argv[++i], "r"))) {   // impurity-profile-spectrum file name
        printf("\nERROR: cannot open impurity profile spectrum file %s\n\n", argv[i]);
        return 1;
      }
      if (strstr(argv[i], ".spe")) fread(setup.rho_z_spe, 36, 1, fp); // skip over .spe-format header
      for (j=0; j<1024; j++) setup.rho_z_spe[j] = 0;
      fread(setup.rho_z_spe, sizeof(setup.rho_z_spe), 1, fp);
      fclose(fp);
    } else if (strstr(argv[i], "-z")) {
      j = atoi(argv[++i]);
      if (j >= 199) {
        printf("Error; rho_spectrum_offset_mm = %d is too large! max = 200\n", j);
        return 1;
      }
      for (k=j; k <200; k++) setup.rho_z_spe[k-j] = setup.rho_z_spe[k];
    } else if (strstr(argv[i], "-g") || strstr(argv[i], "-o")) { // ignore these
      continue;
    } else {
      printf("Possible options:\n"
	     "      -b bias_volts\n"
	     "      -w {0,1,2} (for WV options)\n"
	     "      -p {0,1}   (for WP options)\n"
             "      -r rho_spectrum_file_name\n"
             "      -z <rho_spectrum_offset_mm>\n");
      return 1;
    }
  }

  if (setup.inner_taper_length > setup.hole_length) {
    setup.inner_taper_length = setup.hole_length;
    printf(" Warning: inner_taper_length limited to hole_length = %.1f\n", setup.hole_length);
  }

  if (setup.verbosity >= CHATTY && setup.rho_z_spe[0] != 0) {
    printf(" z(mm)   rho\n");
    for (j=0; j < (int) setup.xtal_length; j++)
      printf(" %3d  %7.3f\n", j, setup.rho_z_spe[j]);
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

  if (setup.write_WP) {
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
  }

  /* -------------- estimate depletion voltage */
  double min = BV, min2 = BV, dV, dW, testv;
  int    vminr=0, vminz=0;
  int    dz[4] = {1, -1, 0, 0}, dr[4] = {0, 0, 1, -1};
  if (setup.write_WP) {
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
        min = min2;
      } else {
        printf("Estimated depletion voltage = %.0f V\n", BV - min);
      }
    }

    printf("Minimum bulk field = %.2f V/cm at (r,z) = (%.1f, %.1f) mm\n\n",
           setup.Emin, setup.rmin, setup.zmin);
  }

#ifdef EMIN_OVERBIAS
  // calculate a new minimum electric field at a specified bias above depletion (usually 500 V)

  if (min > 0) {
    double db = EMIN_OVERBIAS - min;
    int k = 3.0/grid;
    float E, ez, er, r, z;
    // adjust potential
    for (j = 1; j < R; j++) {
      for (i = 1; i < L; i++) {
        setup.vsave[i][j] -= db * setup.v[1][i][j];
      }
    }
    setup.Emin = 99999.9;
    setup.rmin = setup.zmin = 999.9;
    /* calculate new field */
    for (j = 1; j < R-k-1; j++) {
      r = (j-1) * grid;
      for (i = k+1; i < L-k-1; i++) {
        z = (i-1) * grid;
        if (setup.point_type[i][j]     != INSIDE ||
            setup.point_type[i + k][j] != INSIDE ||  // point is at least 3 mm from a boundary
            setup.point_type[i - k][j] != INSIDE ||
            setup.point_type[i][j + k] != INSIDE ||
            (j > k+1 && setup.point_type[i][j - k] != INSIDE)) continue;
        // calc E
        er = 0;
        if (j>1) er = (setup.vsave[i][j-1] - setup.vsave[i][j+1])/(0.2*grid);
        ez = (setup.vsave[i-1][j] - setup.vsave[i+1][j])/(0.2*grid);
        E = sqrt(er*er + ez*ez);
        /* check for minimum field inside bulk of detector */
        if (E > 0.1 && E < setup.Emin) {
          setup.Emin = E;
          setup.rmin = r;
          setup.zmin = z;
        }
      }
    }
    printf("Minimum bulk field at %.0fV (%dV overbias) is %.2f V/cm at (r,z) = (%.1f, %.1f) mm\n\n",
           BV-min+EMIN_OVERBIAS, EMIN_OVERBIAS, setup.Emin, setup.rmin, setup.zmin);
  }
#endif
 
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
  //int    L  = setup->xtal_length/grid + 2.99;
  //int    R  = setup->xtal_radius/grid + 2.99;
  int    L  = lrint(setup->xtal_length/grid)+2;
  int    R  = lrint(setup->xtal_radius/grid)+2;
  float  r, z, E;
  FILE   *file;
  double ***v = setup->v;
  cyl_pt **e;

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
  
  if ((e = (cyl_pt **) malloc(R*sizeof(*e))) == NULL) {
    printf("ERROR: Malloc failed\n");
    fclose(file);
    return 1;
  }
  for (i = 0; i < R; i++){
    if ((e[i] = (cyl_pt *) calloc(L, sizeof(*e[i]))) == NULL) {
      printf("ERROR: Calloc failed\n");
      fclose(file);
      return 1;
    }
  }

  for (j = 1; j < R; j++) {
    r = (j-1) * grid;
    for (i = 1; i < L; i++) {
      z = (i-1) * grid;
      // calc E in r-direction
      if (j == 1) {  // r = 0; symmetry implies E_r = 0
        e[j][i].r = 0;
      } else if (setup->point_type[i][j] == CONTACT_EDGE) {
        e[j][i].r = ((v[new][i][j] - v[new][i][j+1])*setup->dr[1][i][j] +
               (v[new][i][j-1] - v[new][i][j])*setup->dr[0][i][j]) / (0.2*grid);
      } else if (setup->point_type[i][j] < INSIDE &&
                 setup->point_type[i][j-1] == CONTACT_EDGE) {
        e[j][i].r =  (v[new][i][j-1] - v[new][i][j]) * setup->dr[1][i][j-1] / ( 0.1*grid) ;
      } else if (setup->point_type[i][j] < INSIDE &&
                 setup->point_type[i][j+1] == CONTACT_EDGE) {
        e[j][i].r =  (v[new][i][j] - v[new][i][j+1]) * setup->dr[0][i][j+1] / ( 0.1*grid) ;
      } else if (j == R-1) {
        e[j][i].r = (v[new][i][j-1] - v[new][i][j])/(0.1*grid);
      } else {
        e[j][i].r = (v[new][i][j-1] - v[new][i][j+1])/(0.2*grid);
      }
      // calc E in z-direction
      if (setup->point_type[i][j] == CONTACT_EDGE) {
        e[j][i].z = ((v[new][i][j] - v[new][i+1][j])*setup->dz[1][i][j] +
               (v[new][i-1][j] - v[new][i][j])*setup->dz[0][i][j]) / (0.2*grid);
      } else if (setup->point_type[i][j] < INSIDE &&
                 setup->point_type[i-1][j] == CONTACT_EDGE) {
        e[j][i].z =  (v[new][i-1][j] - v[new][i][j]) * setup->dz[1][i-1][j] / ( 0.1*grid) ;
      } else if (setup->point_type[i][j] < INSIDE &&
                 setup->point_type[i+1][j] == CONTACT_EDGE) {
        e[j][i].z =  (v[new][i][j] - v[new][i+1][j]) * setup->dz[0][i+1][j] / ( 0.1*grid) ;
      } else if (i == 1) {
        e[j][i].z = (v[new][i][j] - v[new][i+1][j])/(0.1*grid);
      } else if (i == L-1) {
        e[j][i].z = (v[new][i-1][j] - v[new][i][j])/(0.1*grid);
      } else {
        e[j][i].z = (v[new][i-1][j] - v[new][i+1][j])/(0.2*grid);
      }

      /* temporarily store E in e[j][i].phi */
      E = e[j][i].phi = sqrt(e[j][i].r*e[j][i].r + e[j][i].z*e[j][i].z);
      /* check for minimum field inside bulk of detector */
      int k = 3.0/grid;
      if (E > 0.1 && E < setup->Emin &&
          i > k+1 && j < R-k-1 && i < L-k-1 &&
          setup->point_type[i][j] == INSIDE &&
          setup->point_type[i + k][j] == INSIDE &&  // point is at least 3 mm from a boundary
          setup->point_type[i - k][j] == INSIDE &&
          setup->point_type[i][j + k] == INSIDE &&
          (j < k+1 || setup->point_type[i][j - k] == INSIDE)) {
        setup->Emin = E;
        setup->rmin = r;
        setup->zmin = z;
      }
    }
  }

  if (strstr(setup->field_name, "unf")) {
    fprintf(file, "#\n## start of unformatted data\n");
    i = R-1; j = L-1;
    fwrite(&i, sizeof(int), 1, file);
    fwrite(&j, sizeof(int), 1, file);
    for (i = 1; i < R; i++) {
      for (j = 1; j < L; j++) e[i][j].phi = 0;
      if (fwrite(&e[i][1], sizeof(cyl_pt), L-1, file) != L-1) {
        printf("ERROR while writing %s\n", setup->field_name);
      }
    }
  } else {
    fprintf(file, "#\n## r (mm), z (mm), V (V),  E (V/cm), E_r (V/cm), E_z (V/cm)\n");
    for (j = 1; j < R; j++) {
      r = (j-1) * grid;
      for (i = 1; i < L; i++) {
        z = (i-1) * grid;
        E = e[j][i].phi;
        fprintf(file, "%7.2f %7.2f %7.1f %7.1f %7.1f %7.1f\n",
                r, z, v[new][i][j], E, e[j][i].r, e[j][i].z);
      }
      fprintf(file, "\n");
    }
  }
  fclose(file);
  for (i = 0; i < R; i++) free(e[i]);
  free(e);

  if (!setup->write_WP)
    printf("\n Minimum bulk field = %.2f V/cm at (r,z) = (%.1f, %.1f) mm\n\n",
           setup->Emin, setup->rmin, setup->zmin);

  if (0) { /* write point_type to output file */
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
  //int    L  = setup->xtal_length/grid + 2.99;
  //int    R  = setup->xtal_radius/grid + 2.99;
  int    L  = lrint(setup->xtal_length/grid)+2;
  int    R  = lrint(setup->xtal_radius/grid)+2;
  float  r, z, w;
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

  if (strstr(setup->field_name, "unf")) {
    fprintf(file, "#\n## start of unformatted data\n");
    i = R-1; j = L-1;
    fwrite(&i, sizeof(int), 1, file);
    fwrite(&j, sizeof(int), 1, file);
    for (i = 1; i < R; i++) {
      for (j = 1; j < L; j++) {
        w = setup->v[new][j][i];
        fwrite(&w, sizeof(float), 1, file);
      }
    }
  } else {
    fprintf(file, "#\n## r (mm), z (mm), WP\n");
    for (j = 1; j < R; j++) {
      r = (j-1) * grid;
      for (i = 1; i < L; i++) {
        z = (i-1) * grid;
        fprintf(file, "%7.2f %7.2f %12.6e\n", r, z, setup->v[new][i][j]);
      }
      fprintf(file, "\n");
    }
  }
  fclose(file);

  return 0;
 } /* write_wp */

/* -------------------------------------- dist_from_contact ------------------- */
float dist_from_contact(cyl_pt pt, cyl_pt delta, MJD_Siggen_Setup *setup) {
  float  factor = 1, d = 0.5;
  cyl_pt test;
  int    n;

  for (n=0; n<7; n++) {  // 7 steps => 1/128 precision
    test.r = pt.r + factor * delta.r;
    test.z = pt.z + factor * delta.z;
    if (outside_detector_cyl(test, setup)) {
      factor -= d;
    } else {
      if (n == 0) return -1.0;
      factor += d;
    } 
    d /= 2.0;
  }
  return factor;
} /* dist_from_contact */

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
      (setup->dr[0]      = malloc(L * sizeof(*setup->dr[0])))      == NULL ||
      (setup->dr[1]      = malloc(L * sizeof(*setup->dr[1])))      == NULL ||
      (setup->dz[0]      = malloc(L * sizeof(*setup->dz[0])))      == NULL ||
      (setup->dz[1]      = malloc(L * sizeof(*setup->dz[1])))      == NULL ||
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
        (setup->dr[0][i]      = malloc(R * sizeof(**setup->dr[0])))      == NULL ||
	(setup->dr[1][i]      = malloc(R * sizeof(**setup->dr[1])))      == NULL ||
        (setup->dz[0][i]      = malloc(R * sizeof(**setup->dz[0])))      == NULL ||
	(setup->dz[1][i]      = malloc(R * sizeof(**setup->dz[1])))      == NULL ||
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

  for (i = 1; i < L; i++) {
    for (j = 0; j < R; j++) {
      setup->dr[0][i][j] = setup->s2[j];   //  for r-1
      setup->dr[1][i][j] = setup->s1[j];   //  for r+1
      setup->dz[0][i][j] = 1;       //  for z-1
      setup->dz[1][i][j] = 1;       //  for z+1
    }
  }

  /* set up pixel point types for boundary conditions etc */
  float  d, g = 0.05 * grid, lith  = setup->Li_thickness;
  cyl_pt pt, pt1, pt2, pt3, pt4;

  setup->rmax = setup->xtal_radius - lith;
  setup->zmax = setup->xtal_length - lith;
  setup->hole_radius += lith;
  setup->hole_bullet_radius += lith;
  setup->bottom_taper_length += 0.71*lith; // add top tapers?

  for (i = 1; i < L; i++) {
    pt.z = pt3.z = pt4.z = (i-1) * grid;
    pt1.z = pt.z + g;
    pt2.z = pt.z - g;
    for (j = 1; j < R; j++) {
      pt.r = pt1.r = pt2.r = (j-1) * grid;
      pt3.r = pt.r + g;
      pt4.r = pt.r - g;
      setup->point_type[i][j] = INSIDE;

      /* see if pixel is ouside (or very nearly outside) the detector bulk */
      if (i == 1 || (pt.r >= setup->wrap_around_radius && pt.z - g < lith) ||
          outside_detector_cyl(pt1, setup) || outside_detector_cyl(pt2, setup) ||
          outside_detector_cyl(pt3, setup) || outside_detector_cyl(pt4, setup)) {

        /* check for inside ditch
           boundary condition at Ge-vacuum interface:
           epsilon0 * E_vac = espilon_Ge * E_Ge  */
        if (setup->ditch_depth > 0 && pt.z < setup->ditch_depth + grid &&
            pt.r <= setup->wrap_around_radius &&
            pt.r >= setup->wrap_around_radius - setup->ditch_thickness) {
          setup->point_type[i][j] = DITCH;
          setup->eps[i][j] = setup->eps_dz[i][j] = setup->eps_dr[i][j] = 1.0;

          /* check for inside (point) contact */
        } else if (pt.z < setup->pc_length + g && pt.r < setup->pc_radius+g) {
          setup->point_type[i][j] = PC;

        /* check for passivated area  */
        } else if (i == 1 && pt.r < setup->wrap_around_radius &&      // BEGE/ICPC
                   pt.r < setup->rmax - setup->bottom_taper_length) { // PPC
          setup->point_type[i][j] = PASSIVE;

        /* only remaining surface is HV contact */
        } else  {
          setup->point_type[i][j] = HVC;
        }
      }
    }
  }

  /* find the pixels next to the contact surfaces */
  cyl_pt dp1, dp2, dp3, dp4;
  dp1.z = dp3.r = grid;
  dp2.z = dp4.r = -grid;
  dp1.r = dp2.r = dp3.z = dp4.z = 0;
  for (i = 2; i < L-1; i++) {
    pt.z = (i-1) * grid;
    for (j = 1; j < R; j++) {
      pt.r = (j-1) * grid;
      if (setup->point_type[i][j] == INSIDE &&
          (setup->point_type[i+1][j] < INSIDE || setup->point_type[i-1][j] < INSIDE ||
           setup->point_type[i][j+1] < INSIDE || (j > 1 && setup->point_type[i][j-1] < INSIDE))) {
        setup->point_type[i][j] = CONTACT_EDGE;
        /* find distance to contact surface */
        if (setup->point_type[i+1][j] < INSIDE && (d = dist_from_contact(pt, dp1, setup)) > 0)
          setup->dz[1][i][j] = 1.0/d;
        if (setup->point_type[i-1][j] < INSIDE && (d = dist_from_contact(pt, dp2, setup)) > 0)
          setup->dz[0][i][j] = 1.0/d;
        else if (setup->point_type[i-1][j] < INSIDE &&
                 pt.r >= setup->wrap_around_radius && pt.z - grid < lith)
          setup->dz[0][i][j] = grid/(pt.z - lith);
        if (setup->point_type[i][j+1] < INSIDE && (d = dist_from_contact(pt, dp3, setup)) > 0)
          setup->dr[1][i][j] = setup->s1[j] * 1.0/d;
        if (j > 1 && setup->point_type[i][j-1] < INSIDE &&
            (d = dist_from_contact(pt, dp4, setup)) > 0)
          setup->dr[0][i][j] = setup->s2[j] * 1.0/d;
      }
    }
  }
  setup->hole_radius -= lith;
  setup->hole_bullet_radius -= lith;
  setup->bottom_taper_length -= 0.71*lith;

  /* for pixels adjacent to the ditch, set point_type to DITCH_EDGE
     and for z=0, set flag for passivated surface */
  for (i = 1; i < L; i++) {
    for (j = 1; j < R; j++) {
      setup->eps_dr[i][j-1] = (setup->eps[i][j-1] + setup->eps[i][j]) / 2.0f;
      setup->eps_dz[i-1][j] = (setup->eps[i-1][j] + setup->eps[i][j]) / 2.0f;
      if (setup->point_type[i][j] == INSIDE &&
          (setup->point_type[i-1][j] == DITCH ||
           setup->point_type[i][j-1] == DITCH ||
           setup->point_type[i][j+1] == DITCH)) setup->point_type[i][j] = DITCH_EDGE;
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
  } else if (setup->rho_z_spe[(int) setup->xtal_length] == 0) {
    printf("Error: Impurity profile spectrum is too short? rho[xtal_len = %d] = 0\n", (int) setup->xtal_length);
    return -1;
  } else {
    for (i = 1; i < L; i++)
      imp_z[i] = e_over_E * grid*grid / 4.0 * setup->rho_z_spe[(int) ((i-1) * grid)];
  }
  for (j = 1; j < R; j++) {
    r = (j-1) * grid;
    if (setup->impurity_rpower > 0.1) {
      imp_ra = setup->impurity_radial_add * e_over_E * grid*grid / 4.0 *
        pow((double) r / setup->xtal_radius, setup->impurity_rpower);
      imp_rm = 1.0 + (setup->impurity_radial_mult - 1.0f) *
        pow((double) r / setup->xtal_radius, setup->impurity_rpower);
    }
    for (i = 1; i < L; i++)  setup->impurity[i][j] = imp_z[i] * imp_rm + imp_ra;
    if (setup->point_type[1][j] == PASSIVE) {
      setup->impurity[1][j] += setup->impurity_surface * e_over_E * grid/4.0;
    }
    /* reduce charge volume for CONTACT_EDGE pixels */
    for (i = 1; i < L; i++) {
      if (setup->point_type[i][j] == CONTACT_EDGE) {
        setup->impurity[i][j] /=
          SQ(setup->dz[1][i][j]*setup->dz[0][i][j] * setup->dr[1][i][j]*setup->dr[0][i][j]);
      }
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
  int    L  = lrint(setup->xtal_length/grid)+2;
  int    R  = lrint(setup->xtal_radius/grid)+2;
  double eps_sum, v_sum, dif, sum_dif=0, max_dif;
  double ***v = setup->v, **eps_dr = setup->eps_dr, **eps_dz = setup->eps_dz;
  double ***dr = setup->dr, ***dz = setup->dz;
  double *s1 = setup->s1, *s2 = setup->s2;
  double e_over_E = 11.310; // e/epsilon; for 1 mm2, charge units 1e10 e/cm3, espilon = 16*epsilon0


  if (ev_calc) {
    // for field calculation, save impurity value along passivated surface
    for (r = 1; r < R; r++)
      setup->impurity[0][r] = setup->impurity[1][r];
  } else {
    // for WP calculation, clear all impurity values
    for (z = 0; z < L; z++) {
      for (r = 1; r < R; r++) {
        setup->impurity[z][r] = 0;
      }
    }
  }

  for (iter = 0; iter < setup->max_iterations; iter++) {

    /* the following definition of the factor for over-relaxation improves convergence
           time by a factor ~ 70-120 for a 2kg ICPC detector, grid = 0.1 mm
         OR_fact increases with increasing volxel count (L*R)
               and with increasing iteration number
         0.997 is maximum asymptote for very large pixel count and iteration number */
    double OR_fact;
    if (ev_calc)  OR_fact = (1.991 - 1500.0/(L*R));
    else          OR_fact = (1.992 - 1500.0/(L*R));
    if (OR_fact < 1.4) OR_fact = 1.4;
    // if (iter == 0) printf("OR_fact = %f\n", OR_fact);
    if (iter < 1) OR_fact = 1.0;

    old = new;
    new = 1 - new;
    sum_dif = 0;
    max_dif = 0;

    if (setup->vacuum_gap > 0) {   // modify impurity value along passivated surface
      for (r = 1; r < R; r++)      //   due to surface charge induced by capacitance
        setup->impurity[1][r] = setup->impurity[0][r] -
          v[old][1][r] * 5.52e-4 * e_over_E * grid / setup->vacuum_gap;
    }

    /* start from z=1 and r=1 so that (z,r)=0 can be
       used for reflection symmetry around r=0 or z=0 */
    for (z = 1; z < L; z++) {
      /* manage r=0 reflection symmetry */
      setup->v[old][z][0] = setup->v[new][z][0] = setup->v[old][z][2];

      for (r = 1; r < R; r++) {
        if (setup->point_type[z][r] < INSIDE) continue;   // HV or point contact

        if (setup->point_type[z][r] < DITCH) {       // normal bulk or passivated surface, no complications
          v_sum = (v[old][z+1][r] + v[old][z][r+1]*s1[r] +
                   v[new][z-1][r] + v[new][z][r-1]*s2[r]);
          if (r > 1) eps_sum = 4;
          else       eps_sum = 2 + s1[r] + s2[r];
        } else if (setup->point_type[z][r] == CONTACT_EDGE) {  // adjacent to the contact
          v_sum = (v[old][z+1][r]*dz[1][z][r] + v[old][z][r+1]*dr[1][z][r] +
                   v[new][z-1][r]*dz[0][z][r] + v[new][z][r-1]*dr[0][z][r]);
          eps_sum = dz[1][z][r] + dr[1][z][r] + dz[0][z][r] + dr[0][z][r];
        } else if (setup->point_type[z][r] >= DITCH) {  // in or adjacent to the ditch
          v_sum = (v[old][z+1][r]*eps_dz[z  ][r] + v[old][z][r+1]*eps_dr[z][r  ]*s1[r] +
                   v[new][z-1][r]*eps_dz[z-1][r] + v[new][z][r-1]*eps_dr[z][r-1]*s2[r]);
          eps_sum = (eps_dz[z][r]   + eps_dr[z][r]  *s1[r] +
                     eps_dz[z-1][r] + eps_dr[z][r-1]*s2[r]);
        }

        // calculate the interpolated mean potential and the effect of the space charge

        if ((ev_calc || (setup->vacuum_gap > 0 && z == 1)) &&
            setup->point_type[z][r] < CONTACT_EDGE && r > 1 && z > 1) {   // normal bulk, no complications
          v[new][z][r] = (1.0-OR_fact)*v[old][z][r] + OR_fact * (v_sum / eps_sum + setup->impurity[z][r]);
        } else if (ev_calc || (setup->vacuum_gap > 0 && z == 1)) {
          v[new][z][r] = v_sum / eps_sum + setup->impurity[z][r];
        } else if (setup->point_type[z][r] < CONTACT_EDGE && r > 1 && z > 1) {   // normal bulk, no complications
          v[new][z][r] = (1.0-OR_fact)*v[old][z][r] + OR_fact * v_sum / eps_sum;
        } else {                          // over-relaxation at the edges seems to make things worse
          v[new][z][r] = v_sum / eps_sum;
        }

        // calculate difference from last iteration, for convergence check
        dif = fabs(v[old][z][r] - v[new][z][r]);
        sum_dif += dif;
        if (max_dif < dif) max_dif = dif;
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
    if ( ev_calc && max_dif < 0.0008) break;  // comment out if you want convergence at the numerical error level
    if (!ev_calc && max_dif < 0.0000000001) break;
    if (!ev_calc && max_dif < 0.000001) break;  // comment out if you want convergence at the numerical error level

    /* every 100 iterations, check that detector is really depleted*/
    if (ev_calc && iter > 190 && iter%100 == 0) {
      for (z = 1; z < L; z++) {
        setup->v[old][z][0] = setup->v[old][z][2];
        for (r = 1; r < R; r++) {
          if (setup->point_type[z][r] < INSIDE) continue;   // HV or point contact
          if (v[new][z][r] < 0 ||
              (v[new][z][r] < v[new][z][r] &&
               v[new][z][r] < v[new][z][r] &&
               v[new][z][r] < v[new][z][r] &&
               v[new][z][r] < v[new][z][r])) {
            printf("Detector may not be fully depleted. Switching to ev_relax_undep()\n");
            ev_relax_undep(setup);
            return 0;
          }
        }
      }

    }
  }

  printf(">> %d %.16f\n\n", iter, sum_dif);
  if (setup->vacuum_gap > 0) {   // restore impurity value along passivated surface
    for (r = 1; r < R; r++)
      setup->impurity[1][r] = setup->impurity[0][r];
  }

  return 0;
} /* do_relax */

/* -------------------------------------- ev_relax_undep ------------------- */
/*  This function, unlike do_relax() above, properly handles undepleted detectors.
    Note that this function uses a modified sequential over-relaxtion algorithm,
    while do_relax() above uses a more standard text-book version.
 */
int ev_relax_undep(MJD_Siggen_Setup *setup) {
  int    old = 1, new = 0, iter, r, z, bvn;
  float  grid = setup->xtal_grid;
  int    L  = lrint(setup->xtal_length/grid)+2;
  int    R  = lrint(setup->xtal_radius/grid)+2;
  double eps_sum, v_sum, save_dif, min;
  double dif, sum_dif=0, max_dif, bubble_volts=0;
  double ***v = setup->v, **eps_dr = setup->eps_dr, **eps_dz = setup->eps_dz;
  double ***dr = setup->dr, ***dz = setup->dz;
  double *s1 = setup->s1, *s2 = setup->s2;
  char   **undep = setup->undepleted;
  double e_over_E = 11.310; // e/epsilon; for 1 mm2, charge units 1e10 e/cm3, espilon = 16*epsilon0


  // save impurity value along passivated surface
  for (r = 1; r < R; r++)
    setup->impurity[0][r] = setup->impurity[1][r];

  /* initialise the undepleted array for use with bubble depletion */
  for (z = 1; z < L; z++) {
    for (r = 1; r < R; r++) {
      if (setup->point_type[z][r] >= INSIDE) undep[r][z] = 0;
    }
  }

  for (iter = 0; iter < setup->max_iterations; iter++) {

    double OR_fact = ((0.997 - 300.0/(L*R)) * (1.0 - 0.9/(double)(1+iter/6)));
    if (300.0/(L*R) > 0.5) OR_fact = (0.5 * (1.0 - 0.9/(double)(1+iter/6)));
    if (iter < 2) OR_fact = 0.0;

    old = new;
    new = 1 - new;
    sum_dif = 0;
    max_dif = 0;
    bubble_volts = 0;
    bvn = 0;

    if (setup->vacuum_gap > 0) {   // modify impurity value along passivated surface
      for (r = 1; r < R; r++)      //   due to surface charge induced by capacitance
        setup->impurity[1][r] = setup->impurity[0][r] -
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
          if (r > 1) eps_sum = 4;
          else       eps_sum = 2 + s1[r] + s2[r];
        } else if (setup->point_type[z][r] == CONTACT_EDGE) {  // adjacent to the contact
          v_sum = (v[old][z+1][r]*dz[1][z][r] + v[old][z][r+1]*dr[1][z][r] +
                   v[old][z-1][r]*dz[0][z][r] + v[old][z][r-1]*dr[0][z][r]);
          eps_sum = dz[1][z][r] + dr[1][z][r] + dz[0][z][r] + dr[0][z][r];
        } else if (setup->point_type[z][r] >= DITCH) {  // in or adjacent to the ditch
          v_sum = (v[old][z+1][r]*eps_dz[z  ][r] + v[old][z][r+1]*eps_dr[z][r  ]*s1[r] +
                   v[old][z-1][r]*eps_dz[z-1][r] + v[old][z][r-1]*eps_dr[z][r-1]*s2[r]);
          eps_sum = (eps_dz[z][r]   + eps_dr[z][r]  *s1[r] +
                     eps_dz[z-1][r] + eps_dr[z][r-1]*s2[r]);
        }

        // calculate the interpolated mean potential and the effect of the space charge
        min = fminf(fminf(v[old][z+1][r], v[old][z][r+1]),
                    fminf(v[old][z-1][r], v[old][z][r-1]));
        v[new][z][r] = v_sum / eps_sum + setup->impurity[z][r];

        undep[r][z] /= 2;
        if (v[new][z][r] <= 0) {
          v[new][z][r] = 0;
          undep[r][z] = 4;  // do not do over-relaxation for 3 iterations
        } else if (v[new][z][r] <= min) {
          if (bubble_volts == 0) bubble_volts = min + 0.2*grid*grid; // finer grids require smaller increment here
          v[new][z][r] = bubble_volts;
          bvn++;
          undep[r][z] = 8;  // do not do over-relaxation for 4 iterations
        }

        // calculate difference from last iteration, for convergence check
        dif = v[old][z][r] - v[new][z][r];
        if (dif < 0) dif = -dif;
        sum_dif += dif;
        if (max_dif < dif) max_dif = dif;
        // do over-relaxation
        if (!undep[r][z])  v[new][z][r] += OR_fact*save_dif;
      }
    }

    // report results for some iterations
    if (iter < 10 || (iter < 600 && iter%100 == 0) || iter%1000 == 0) {
      if (0) {
        printf("%5d %d %d %.10f %.10f\n", iter, old, new, max_dif, sum_dif/(L-2)/(R-2));
      } else {
        printf("%5d %d %d %.10f %.10f ; %.10f %.10f bubble %.2f %d\n",
               iter, old, new, max_dif, sum_dif/(L-2)/(R-2),
               v[new][L/2][R/2], v[new][L/3][R/3], bubble_volts, bvn);
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
      if (setup->point_type[z][r] < INSIDE) {
        undep[r][z] = ' ';
      } else if (undep[r][z] == 0) {
        undep[r][z] = '.';
      } else {
        if (undep[r][z] > 4) undep[r][z] = 'B';  // identifies pinch-off
        else undep[r][z] = '*';
        setup->fully_depleted = 0;
      }
    }
  }

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
  int    L  = lrint(setup->xtal_length/grid)+2;
  int    R  = lrint(setup->xtal_radius/grid)+2;
  double eps_sum, v_sum, save_dif, pinched_sum1, pinched_sum2;
  double dif, sum_dif=0, max_dif;
  double ***v = setup->v, **eps_dr = setup->eps_dr, **eps_dz = setup->eps_dz;
  double ***dr = setup->dr, ***dz = setup->dz;
  double *s1 = setup->s1, *s2 = setup->s2;
  double e_over_E = 11.310; // e/epsilon; for 1 mm2, charge units 1e10 e/cm3, espilon = 16*epsilon0


  for (iter = 0; iter < setup->max_iterations; iter++) {

   double OR_fact = ((0.997 - 300.0/(L*R)) * (1.0 - 0.9/(double)(1+iter/6)));
    if (300.0/(L*R) > 0.5) OR_fact = (0.5 * (1.0 - 0.9/(double)(1+iter/6)));
    if (iter < 2) OR_fact = 0.0;

    old = new;
    new = 1 - new;
    sum_dif = 0;
    max_dif = 0;
    pinched_sum1 = pinched_sum2 = 0.0;

    if (setup->vacuum_gap > 0) {   // modify impurity value along passivated surface
      for (r = 1; r < R; r++)      //   due to surface charge induced by capacitance
        setup->impurity[1][r] = -v[old][1][r] * 5.52e-4 * e_over_E * grid / setup->vacuum_gap;
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
          if (r > 1) eps_sum = 4;
          else       eps_sum = 2 + s1[r] + s2[r];
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
        } else if (setup->point_type[z][r] == CONTACT_EDGE) {  // adjacent to the contact
          v_sum = (v[old][z+1][r]*dz[1][z][r] + v[old][z][r+1]*dr[1][z][r] +
                   v[old][z-1][r]*dz[0][z][r] + v[old][z][r-1]*dr[0][z][r]);
          eps_sum = dz[1][z][r] + dr[1][z][r] + dz[0][z][r] + dr[0][z][r];
        } else if (setup->point_type[z][r] >= DITCH) {  // in or adjacent to the ditch
          v_sum = (v[old][z+1][r]*eps_dz[z  ][r] + v[old][z][r+1]*eps_dr[z][r  ]*s1[r] +
                   v[old][z-1][r]*eps_dz[z-1][r] + v[old][z][r-1]*eps_dr[z][r-1]*s2[r]);
          eps_sum = (eps_dz[z][r]   + eps_dr[z][r]  *s1[r] +
                     eps_dz[z-1][r] + eps_dr[z][r-1]*s2[r]);
        }

        if (setup->point_type[z][r] != PINCHOFF) {
          // calculate the interpolated mean potential and the effect of the space charge
          if (setup->vacuum_gap > 0 && z == 1)
            v[new][z][r] = v_sum / eps_sum + setup->impurity[z][r];
          else
            v[new][z][r] = v_sum / eps_sum;

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
      for (z=1; z<L; z++) {
        for (r=1; r<R; r++) {
          if (setup->point_type[z][r] == PINCHOFF) {
            v[new][z][r] = pinched_sum1 / pinched_sum2;
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
