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
      - When refining the adaptive grid size, if the potential is sufficiently linear,
      and the grid size is a small fraction of the crystal size,
      then try to fix the potential to the interpolated value from the previous grid-size step.
      Note that this option reduces the accuracy of the minimum field and depletion voltage.
      - Removed -c flag for config files; config file should always be first argument

   TO DO:
      - add other bulletizations
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

int report_config(FILE *fp_out, char *config_file_name);
int interpolate(MJD_Siggen_Setup *setup, double **v[2], float grid);
int interpolate2(double **v[2], int **bulk, int ig, int L, int R);
int set_boundary_conditions(MJD_Siggen_Setup *setup, float grid, double **v[2],
                            int **bulk, double **vfraction, float BV);
int do_relax(int R, int L, int LC, int istep, int max_its, double **v[2], int **bulk, int ig,
             double **eps_dr, double **eps_dz, double *s1, double *s2,
             char **undepleted, float **imp, int evcalc, int **d[4]);
int outside_detector(float z, float r, MJD_Siggen_Setup *setup);



int main(int argc, char **argv)
{

  MJD_Siggen_Setup setup;

  /* --- default values, normally over-ridden by values in a *.conf file --- */
  int   R  = 0;  // radius of detector, in grid lengths
  int   L  = 0;  // length of detector, in grid lengths
  int   RC = 0;  // radius of central contact, in grid lengths
  int   LC = 0;  // length of central contact, in grid lengths
  int   RO = 0;  // radius of wrap-around outer contact, in grid lengths
  int   LO = 0;  // length of ditch next to wrap-around outer contact, in grid lengths
  int   WO = 0;  // width of ditch next to wrap-around outer contact, in grid lengths

  float BV = 0;  // bias voltage
  float N = 1;   // charge density at z=0 in units of e+10/cm3
  float M = 0;   // charge density gradient, in units of e+10/cm4

  int   WV = 0;  // 0: do not write the V and E values to ppc_ev.dat
                 // 1: write the V and E values to ppc_ev.dat
                 // 2: write the V and E values for both +r, -r (for gnuplot, NOT for siggen)
  int   WP = 0;  // 0: do not calculate the weighting potential
                 // 1: calculate the WP and write the values to ppc_wp.dat
  int   WD = 0;  // 0: do not write out depletion surface
                 // 1: write out depletion surface to depl_<HV>.dat
  /* ---  ---  ---  ---  ---  ---  ---  ---  ---  ---  ---  ---  ---  --- */

  double **v[2], **eps, **eps_dr, **eps_dz, **vfraction, *s1, *s2, **vsave;
  char   **undepleted, config_file_name[256], fname[256];
  int    **bulk, ***d;
  float  **imp;

  double e_over_E = 11.31; // e/epsilon; for 1 mm2, charge units 1e10 e/cm3, espilon = 16*epsilon0
  double esum, esum2, pi=3.14159, Epsilon=(8.85*16.0/1000.0);  // permittivity of Ge in pF/mm
  double *imp_ra, *imp_rm, *imp_z, S=0;
  double min, min2, dVn[4], dWn[4], test;
  float  rho_z[256] = {0};
  float  a, b, c, grid = 0.5, grid0;
  float  E, E_r, E_z, bubble_volts=0, gridstep[10] = {0}, Emin, rmin, zmin;
  int    fully_depleted=0, LL=L, RR=R, vminr=0, vminz=0;
  int    i, j, r, z, old, new=0, istep, max_its, ig;
  int    vminr2=0, vminz2=0, fix_adaptive = 0;
  time_t t0=0, t1, t2=0;
  FILE   *file, *fp;


  if (argc < 2 || argc%2 != 0 || read_config(argv[1], &setup)) {
    printf("Usage: %s <config_file_name> [options]\n"
           "   Possible options:\n"
	   "      -b bias_volts\n"
	   "      -w {0,1}  do_not/do write the field file)\n"
	   "      -d {0,1}  do_not/do write the depletion surface)\n"
	   "      -p {0,1}  do_not/do write the WP file)\n"
           "      -f {0,1}  do_not/do try to fix potential values from larger grid sizes)\n"
           "                - note that this reduces accuracy of min field and depl voltage\n"
           "      -r rho_spectrum_file_name\n", argv[0]);
    return 1;
  }
  strncpy(config_file_name, argv[1], sizeof(config_file_name));

  if (setup.xtal_grid < 0.001) setup.xtal_grid = 0.5;
  grid = grid0 = setup.xtal_grid;

  L  = LL = lrint(setup.xtal_length/grid);
  R  = RR = lrint(setup.xtal_radius/grid);
  LC = lrint(setup.pc_length/grid);
  RC = lrint(setup.pc_radius/grid);
  RO = lrint(setup.wrap_around_radius/grid);
  LO = lrint(setup.ditch_depth/grid);
  WO = lrint(setup.ditch_thickness/grid);

  N  = setup.impurity_z0;
  M  = setup.impurity_gradient;
  BV = setup.xtal_HV;
  WV = setup.write_field;
  WP = setup.write_WP;

  for (i=2; i<argc-1; i++) {
    if (strstr(argv[i], "-b")) {
      BV = atof(argv[++i]);   // bias volts
    } else if (strstr(argv[i], "-w")) {
      WV = atoi(argv[++i]);   // write-out options
    } else if (strstr(argv[i], "-d")) {
      WD = atoi(argv[++i]);   // write-out options
    } else if (strstr(argv[i], "-p")) {
      WP = atoi(argv[++i]);   // weighting-potential options
    } else if (strstr(argv[i], "-f")) {
      fix_adaptive = atoi(argv[++i]);   // adaptive-grid-potential fixing option
    } else if (strstr(argv[i], "-r")) {
      if (!(fp = fopen(argv[++i], "r"))) {   // impurity-profile-spectrum file name
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
  }

  if (L*R > 2500*2500) {
    printf("Error: Crystal size divided by grid size is too large!\n");
    return 1;
  }
  if (WV < 0 || WV > 2) WV = 0;

  if (setup.verbosity >= CHATTY) {
    /* give details of detector geometry */
    /*
    printf("\n\n"
           "      Crystal: Radius x Length: %.1f x %.1f mm\n"
           " Bottom taper: %.1f mm\n",
	   grid * (float) R, grid * (float) L, grid * (float) LT);
    if (LH > 0) {
      if (HTL > 0)
        printf("    Core hole: Radius x length: %.1f x %.1f mm,"
               " taper %.1f x %.1f mm (%2.f degrees)\n",
               grid * (float) RH, grid * (float) LH,
               grid * (float) HTW, grid * (float) HTL, setup.taper_angle);
      else
        printf("    Core hole: Radius x length: %.1f x %.1f mm\n",
               grid * (float) RH, grid * (float) LH);
    }
    if (OTL > 0) {
      printf("Outside taper: %.1f mm over %.1f mm (%.2f degrees)\n\n",
             grid * (float) OTW, grid * (float) OTL, setup.taper_angle);
    } else if (LH <= 0) {
      printf("  No core hole or outside taper.\n");
    } else {
      printf("  No outside taper.\n");
    }
    if (LiT > 0) {
      printf("Li contact thickness: %.1f mm\n", grid * (float) LiT);
    }
    if (LC > 0 && setup.bulletize_PC) {
      printf("Point contact: Radius x length: %.1f x %.1f mm, bulletized\n",
             grid * (float) RC, grid * (float) LC);
    } else if (LC > 0) {
      printf("Point cContact: Radius x length: %.1f x %.1f mm, not bulletized\n",
             grid * (float) RC, grid * (float) LC);
    } else {
      printf("Point contact: Radius x length: %.1f x %.1f mm\n",
             grid * (float) RC, grid * (float) LC);
    }
    if (RO <= 0.0 || RO >= R) {
      RO = R - LT;    // inner radius of bottom taper, in grid lengths
      printf(" No wrap-around contact or ditch...\n");
    } else {
      printf("  Wrap-around: Radius x ditch x gap:  %.1f x %.1f x %.1f mm\n",
             grid * (float) RO, grid * (float) LO, grid * (float) WO);
    }
    printf("         Bias: %.0f V\n"
           "   Impurities: (%.3f + %.3fz) e10/cm3\n\n",
           BV, N, M);
    */
  }
    
  if ((BV < 0 && N < 0) || (BV > 0 && N > 0)) {
    printf("ERROR: Expect bias and impurity to be opposite sign!\n");
    return 1;
  } 
  if (N > 0) {
    // swap polarity for n-type material; this lets me assume all voltages are positive
    BV = -BV;
    M = -M;
    N = -N;
  }

  /*
    If grid is too small compared to the crystal size, then it will take too
    long for the relaxation to converge. In that case, we use an adaptive
    grid, where we start out coarse and then refine the grid.
  */
  /*
    DCR Nov 2019 : change to using four simple factors of 2 in grid size,
    rather than earlier more complicated system with fewer but maybe larger steps.
  */
  if (0) {
    gridstep[0] = grid;
    gridstep[1] = 0;
  } else {
    gridstep[0] = grid*9.0;
    gridstep[1] = grid*3.0;
    gridstep[2] = grid;
    gridstep[3] = 0;
  }
  printf("\n---------------- Three grid sizes: %.4f %.4f %.4f\n",
         gridstep[0], gridstep[1], gridstep[2]);
  
  /* ------------------------------------------------------------ */
  /* malloc arrays
     float v[2][L+5][R+5];
     float eps[L+1][R+1], eps_dr[L+1][R+1], eps_dz[L+1][R+1];
     float vfraction[L+1][R+1], s1[R], s2[R], drrc[LC+2], drrc[LC+2], drrh[L+1], drrh[L+1];
     char  undepleted[R+1][L+1];
     int   bulk[L+1][R+1], rrc[LC+2], rrh[L+1];
  */
  // here ig is the number of grid-point distances for the coarsest grid we will use
  ig = lrint(gridstep[0]/grid);

#define ERR { printf("Malloc failed; j = %d\n", j); return 1; }
  printf(" ----------------- ig = %d -----------------\n", ig);
  j = -1; // for ERRor reporting
  if ((v[0]   = malloc((L+ig)*sizeof(*v[0]))) == NULL ||
      (v[1]   = malloc((L+ig)*sizeof(*v[1]))) == NULL ||
      (d      = malloc((L+ig)*sizeof(*d)))    == NULL ||
      (imp    = malloc((L+ig)*sizeof(*imp)))  == NULL ||
      (eps    = malloc((L+ig)*sizeof(*eps)))  == NULL ||
      (eps_dr = malloc((L+ig)*sizeof(*eps_dr))) == NULL ||
      (eps_dz = malloc((L+ig)*sizeof(*eps_dz))) == NULL ||
      (bulk   = malloc((L+ig)*sizeof(*bulk)))   == NULL ||
      (vfraction  = malloc((L+ig)*sizeof(*vfraction)))  == NULL ||
      (undepleted = malloc((R+ig)*sizeof(*undepleted))) == NULL ||
      (imp_ra = malloc((R+ig)*sizeof(*imp_ra))) == NULL ||
      (imp_rm = malloc((R+ig)*sizeof(*imp_rm))) == NULL ||
      (imp_z = malloc((L+ig)*sizeof(*imp_z))) == NULL ||
      (s1 = malloc((R+ig)*sizeof(*s1))) == NULL ||
      (s2 = malloc((R+ig)*sizeof(*s2))) == NULL ||
      (vsave = malloc((L)*sizeof(*vsave))) == NULL) ERR;
  for (j=0; j<L+ig; j++) {
    if ((v[0][j] = malloc((R+ig)*sizeof(**v[0]))) == NULL ||
        (v[1][j] = malloc((R+ig)*sizeof(**v[1]))) == NULL ||
        (d[j]    = malloc((R+ig)*sizeof(**d)))    == NULL ||
        (imp[j]  = malloc((R+ig)*sizeof(**imp)))  == NULL ||
        (eps[j]  = malloc((R+ig)*sizeof(**eps)))  == NULL ||
        (eps_dr[j] = malloc((R+ig)*sizeof(**eps_dr))) == NULL ||
        (eps_dz[j] = malloc((R+ig)*sizeof(**eps_dz))) == NULL ||
        (bulk[j] = malloc((R+ig)*sizeof(**bulk))) == NULL ||
        (vfraction[j] = malloc((R+ig)*sizeof(**vfraction))) == NULL) ERR;
    for (i=0; i<R+ig; i++)
      if ((d[j][i] = malloc(4*sizeof(***d))) == NULL) ERR;
  }
  for (j=0; j<L; j++)
    if ((vsave[j]  = malloc((R)*sizeof(**vsave)))  == NULL) ERR;
  for (j=0; j<R+1; j++) {
    if ((undepleted[j] = malloc((L+1)*sizeof(**undepleted))) == NULL) ERR;
    memset(undepleted[j], ' ', (L+1)*sizeof(**undepleted));
  }
  for (r=0; r<R+ig; r++) {
    imp_ra[r] = 0.0;
    imp_rm[r] = 1.0;
  }

  /* to be safe, initialize overall potential to bias voltage */
  for (z=0; z<L+ig; z++) {
    for (r=0; r<R+ig; r++) {
      v[0][z][r] = v[1][z][r] = BV;
    }
  }
  if (setup.verbosity >= CHATTY) t0 = t2 = time(NULL);  // for calculating elapsed time later
  max_its = MAX_ITS;
  if (setup.max_iterations > 0) max_its = setup.max_iterations;

  /* ----------------------------------------------------------------------------
     now set up and perform the EV relaxation for each of the grid step sizes in turn
  */
  for (istep = 0; istep < 4 && gridstep[istep] > 0; istep++) {
    if (istep > 0) {
      /* not the first go-around, so the previous calculation was on a coarser grid...
	 now interpolate the potential to the final fine grid */
      interpolate(&setup, v, grid);
    }
    grid = gridstep[istep]; // grid size for this go-around
    old = 1;
    new = 0;
    ig = lrintf(grid / grid0);
    /*  e/espilon * area of pixel in mm2 / 4
	for 1 mm2, charge units 1e10 e/cm3, espilon = 16*epsilon0
	4.0 = surface area / volume of voxel in cylindrical (2D) coords
	(this would be 6.0 in cartesian coords) */
    e_over_E = 11.31 * grid*grid / 4.0;
    S = setup.impurity_surface * e_over_E / grid;
    if (setup.verbosity >= NORMAL) printf("grid = %f  ig = %d\n", grid, ig);

    // set up impurity concentration
    if (rho_z[0] != 0) {
      for (z=0; z<L+ig; z++) {
        imp_z[z] = rho_z[(int)(grid0 * (double) z + 0.5)] * e_over_E;
      }
    } else {
      for (z=0; z<L+ig; z++) {
        imp_z[z] = (N + 0.1 * M * grid0 * (double) z + 
                    setup.impurity_quadratic * (1.0 - (double) ((z-L/2)*(z-L/2)) /
                                                (double) (L*L/4))) * e_over_E;
      }
    }
    if (setup.impurity_rpower > 0.1) {
      for (r=0; r<R+ig; r++) {
	imp_ra[r] = setup.impurity_radial_add * e_over_E *
	  pow((double) r / (double) R, setup.impurity_rpower);
	imp_rm[r] = 1.0 + (setup.impurity_radial_mult - 1.0f) *
	  pow((double) r / (double) R, setup.impurity_rpower);
      }
    }
    if (istep == 0) {
      // no previous coarse relaxation, so make initial wild guess at potential:
      for (z=0; z<L; z++) {
	a = BV * (float) (z) / (float) L;
	for (r=0; r<R; r++) {
	  v[0][z][r] = v[1][z][r] = a + (BV - a) * (float) (r) / (float) R;
	}
      }
    }

    /* ------------------------------------------------------------ */
    /* boundary conditions and permittivity
       boundary condition at Ge-vacuum interface:
       epsilon0 * E_vac = espilon_Ge * E_Ge
    */
    set_boundary_conditions(&setup, grid, v, bulk, vfraction, BV);
    /* In the following we divide areas and volumes by pi
       r_bin   rmax  A_top A_outside A_inside  volume  total_surf  out/top  tot/vol
       0     1/2    1/4      1         0       1/4      1.5         4        6  << special case
       1     3/2      2      3         1        2        8        3/2        4
       2     5/2      4      5         3        4       16        5/4        4
       3     7/2      6      7         5        6       24        7/6        4
       r   r+0.5     2r    2r+1      2r-1      2r       8r     (2r+1)/2r     4
       = 1+0.5/r
    */
    // weighting values for the relaxation alg. as a function of r
    s1[0] = 4.0;
    s2[0] = 0.0;
    for (r=ig; r<R+ig; r+=ig) {
      s1[r] = 1.0 + 0.5 / (double) (r / ig);   //  for r+1
      s2[r] = 1.0 - 0.5 / (double) (r / ig);   //  for r-1
    }
    for (z=0; z<L+ig; z+=ig) {
      for (r=0; r<R+ig; r+=ig) {
	eps[z][r] = eps_dz[z][r] = eps_dr[z][r] = 16;   // permittivity inside Ge
	if (z < LO  && r < RO && r > RO-WO-1) eps[z][r] =  1;  // permittivity inside vacuum
	if (z > 0) eps_dz[z-ig][r] = (eps[z-ig][r]+eps[z][r])/2.0f;
	if (r > 0) eps_dr[z][r-ig] = (eps[z][r-ig]+eps[z][r])/2.0f;
      }
    }

  /* ------------------------------------------------------------ */
    // now do the actual relaxation
    for (z=0; z<L; z+=ig) {
      for (r=0; r<R; r+=ig) {
        imp[z][r] = vfraction[z][r] * (imp_z[z]*imp_rm[r] + imp_ra[r]);
        if (r == 0) imp[z][r] /=  1.5;  // special case where volume of voxel is 1/6 of area, not 1/4
        if ((z == 0 && r > RC && r < RO-WO) ||         // passivated surface at z = 0
            (z < LO && (r == RO || r == RO-WO-ig)) ||  // passivated surface on sides of ditch
            (z == LO && r <= RO && r >= RO-WO-ig))     // passivated surface at top of ditch
          imp[z][r] += vfraction[z][r] * S;
        d[z][r][0] = z-ig;
        d[z][r][1] = z+ig;
        d[z][r][2] = r-ig;
        d[z][r][3] = r+ig;
      }
    }
    if (fix_adaptive && istep > 1) {
      // see if the pixel is at least four grid steps away from a boundary (or r=0)
      // remember that ig is the curent pixel size (in units of grid0)
      for (z=42*ig; z<L-4*ig; z+=3*ig) {
        for (r=42*ig; r<R-4*ig; r+=3*ig) {
          //if (r*grid0 < 11) continue;
          if (bulk[z][r] || bulk[z][r-ig] ||
              bulk[z-ig][r] || bulk[z-ig][r-ig] || bulk[z-ig][r+ig]) continue;
          int ok = 1;
          for (i=-4*ig; i<=4*ig; i+=ig) {
            for (j=-4*ig; j<=4*ig; j+=ig) {
              if (bulk[z+i][r+j] &&
                  bulk[z+i][r+j]!= -4 &&
                  bulk[z+i][r+j]!=  4) ok=0; // no, it's not
            }
          }
          if (ok) {   // yes it is
            d[z][r][0] = z-3*ig;
            d[z][r][1] = z+3*ig;
            d[z][r][2] = r-3*ig;
            d[z][r][3] = r+3*ig;
            for (i=-ig; i<=ig; i++) {
              for (j=-ig; j<=ig; j++) {
                bulk[z+i][r+j] = -4;
              }
            }
            if (fix_adaptive > 2) {
              bulk[z][r] = 4;
              imp[z][r] *= 9;
            }
          }
        }
      }
    }
    int n=0;
    for (z=0; z<L; z+=ig) {
      for (r=0; r<R; r+=ig) {
        if (bulk[z][r] == -4 || bulk[z][r] == 4) n++;
      }
    }
    printf("\n ------ n = %d of %d\n", n, L*R/ig/ig);

    do_relax(R, L, LC, istep, max_its, v, bulk, ig, eps_dr, eps_dz,
             s1, s2, undepleted, imp, 1, d);

    fully_depleted = 1;
    for (r=0; r<R+1; r++) {
      for (z=0; z<L+1; z++) {
	if (undepleted[r][z] == '*') {
	  fully_depleted = 0;
	  if (v[new][z][r] > 0.001) undepleted[r][z] = 'B';  // identifies pinch-off
	}
      }
    }
    if (fully_depleted) {
      printf("Detector is fully depleted.\n");
    } else {
      printf("Detector is not fully depleted.\n");
      if (bubble_volts > 0.0f) printf("Pinch-off bubble at %.0f V potential\n", bubble_volts);
    }
    if (setup.verbosity >= CHATTY) {
      t1 = time(NULL);
      printf("\n ^^^^^^^^^^^^^ %d (%d) s elapsed ^^^^^^^^^^^^^^\n",
	     (int) (t1 - t0), (int) (t1 - t2));
      t2 = t1;
    }

    if (istep == 0) {
      // can reduce # of iterations after first go-around
      max_its /= MAX_ITS_FACTOR;
      // report V and E along the axes r=0 and z=0
      if (setup.verbosity >= CHATTY) {
	printf("  z(mm)(r=0)      V   E(V/cm) |  r(mm)(z=0)      V   E(V/cm)\n");
	a = b = v[new][0][0];
	for (z=0; z<L+ig; z+=ig) {
	  printf("%10.1f %8.1f %8.1f  |",
		 ((float) z)*grid0, v[new][z][0], (v[new][z][0] - a)/(0.1*grid));
	  a = v[new][z][0];
	  if (z > R) {
	    printf("\n");
	  } else {
	    r = z;
	    printf("%10.1f %8.1f %8.1f\n",
		   ((float) r)*grid0, v[new][0][r], (v[new][0][r] - b)/(0.1*grid));
	    b = v[new][0][r];
	  }
	}
      }
      // write a little file that shows any undepleted voxels in the crystal
      file = fopen("undepleted.txt", "w");
      for (r=R; r>=0; r--) {
	undepleted[r][L] = '\0';
	fprintf(file, "%s\n", undepleted[r]);
      }
      fclose(file);
    }
  }
  if (grid > grid0) {
    for (z=0; z<L; z+=ig) {
      for (r=0; r<R; r+=ig) {
        for (i=0; i<ig; i++) {
          for (j=0; j<ig; j++) {
            v[new][z+i][r+j] = v[new][z][r];
          }
        }
      }
    }
  }
  if (grid > grid0) interpolate(&setup, v, grid);
  if (fix_adaptive == 2 || fix_adaptive == 4) {
    for (z=0; z<L; z++) {
      for (r=0; r<R; r++) {
        if (bulk[z][r] == -4 || bulk[z][r] == 4) {
          bulk[z][r] = 0;
          d[z][r][0] = z-ig;
          d[z][r][1] = z+ig;
          d[z][r][2] = r-ig;
          d[z][r][3] = r+ig;
          imp[z][r] = vfraction[z][r] * (imp_z[z]*imp_rm[r] + imp_ra[r]);
        }
      }
    }   
    do_relax(R, L, LC, istep, 1000, v, bulk, ig, eps_dr, eps_dz,
             s1, s2, undepleted, imp, 1, d);
  }

  /* ------------------------------------------------------------ */
  if (WV) {
    if (setup.impurity_z0 > 0) {
      // swap voltages back to negative for n-type material
      for (r=0; r<R+1; r++) {
	for (z=0; z<L+1; z++) {
	  v[new][z][r] = -v[new][z][r];
	}
      }
    }
    // write potential and field to output file
    if (!(file = fopen(setup.field_name, "w"))) {
      printf("ERROR: Cannot open file %s for electric field...\n", setup.field_name);
      return 1;
    } else {
      printf("Writing electric field data to file %s\n", setup.field_name);
    }
    /* copy configuration parameters to output file */
    report_config(file, config_file_name);
    fprintf(file, "#\n# HV bias in fieldgen: %.1f V\n", BV);
    if (fully_depleted) {
      fprintf(file, "# Detector is fully depleted.\n");
    } else {
      fprintf(file, "# Detector is not fully depleted.\n");
      if (bubble_volts > 0.0f) fprintf(file, "# Pinch-off bubble at %.0f V potential\n", bubble_volts);
    }
    fprintf(file, "#\n## r (mm), z (mm), V (V),  E (V/cm), E_r (V/cm), E_z (V/cm)\n");

    Emin = 9999.9;
    rmin = zmin = 99.9;
    int RS = 0;
    if (WV > 1) RS = -R;
    for (int rr=RS; rr<R+1; rr++) {
      r = rr;
      if (rr < 0) r= -rr;
      for (z=0; z<L+1; z++) {
	// calc E in r-direction
	if (r==0) {
	  // E_r = (v[new][z][r] - v[new][z][r+1])/(0.1*grid0);
	  E_r = 0;
	} else if (r==R) {
	  E_r = (v[new][z][r-1] - v[new][z][r])/(0.1*grid0);
	} else {
	  E_r = (v[new][z][r-1] - v[new][z][r+1])/(0.2*grid0);
	}
	// calc E in z-direction
	if (z==0) {
	  E_z = (v[new][z][r] - v[new][z+1][r])/(0.1*grid0);
	} else if (z==L) {
	  E_z = (v[new][z-1][r] - v[new][z][r])/(0.1*grid0);
	} else {
	  E_z = (v[new][z-1][r] - v[new][z+1][r])/(0.2*grid0);
	}
        E = sqrt(E_r*E_r + E_z*E_z);
        if (E > 0.1 && E < Emin &&
            (R-r)*grid0 > 5.0 &&                   // more than 5 mm from outer radius
            (L-z)*grid0 > 5.0 && z*grid0 > 5.0 &&  // more than 5 mm from top & bottom
            !outside_detector(grid0 * z - 1.0, grid0 * r, &setup) &&  // outside point contact
            !outside_detector(grid0 * z, grid0 * r - 5.0, &setup) &&  // outside inner well
            !outside_detector(grid0 * z, grid0 * r - 5.0, &setup)) {  // outside any outer taper
          Emin = E;
          rmin = r*grid0;
          zmin = z*grid0;
        }

        fprintf(file, "%7.2f %7.2f %7.1f %7.1f %7.1f %7.1f\n",
                ((float) rr)*grid0,  ((float) z)*grid0, v[new][z][r], E, E_r, E_z);
      }
      fprintf(file, "\n");
    }
    fclose(file);
    printf("\n Minimum bulk field = %.1f V/cm at (r,z) = (%.1f, %.1f) mm\n\n",
           Emin, rmin, zmin);
  }

  /* ------------------------------------------------------------ */
  if (WD) {
    // write depletion surface to output file
    sprintf(fname, "depl_%4.4dV.dat", (int) BV);
    if (!(file = fopen(fname, "w"))) {
      printf("ERROR: Cannot open file %s for depletion...\n", fname);
      return 1;
    } else {
      printf("Writing depletion surface data to file %s\n", fname);
    }
    fprintf(file, "#\n## r (mm), z (mm), \n");
    rmin = zmin = 99.9;
    int RS = -R;
    for (int rr=RS; rr<R+1; rr++) {
      r = rr;
      if (rr < 0) r= -rr;
      for (z=0; z<L+1; z++) {
        j = 0;
        if (undepleted[r][z] == '.') j = 1;
        if (undepleted[r][z] == '*') j = 2;
        if (undepleted[r][z] == 'B') j = 2;
	fprintf(file, "%7.2f %7.2f %d\n",
		((float) rr)*grid,  ((float) z)*grid, j);
      }
      fprintf(file, "\n");
    }
    fclose(file);
  }

  /* ----------------------------------------------------------------------------
     --- weighting potential --- weighting potential --- weighting potential --- 
  */
  if (WP == 0) return 0;
  if (fully_depleted) {
    /* save potential close to point contact,
       to use later when calculating depletion voltage */
    //for (z=0; z<LC+2; z++) {
    //  for (r=0; r<RC+2; r++) {
    for (z=0; z<L; z++) {
      for (r=0; r<R; r++) {
        vsave[z][r] = fabs(v[new][z][r]);
      }
    }
  }

  /*
    -------------------------------------------------------------------------
    now calculate the weighting potential for the central contact
    the WP is also needed for calculating the capacitance
    -------------------------------------------------------------------------
  */

  printf("\nCalculating weighting potential...\n\n");
  if (setup.verbosity >= CHATTY) t0 = t2 = time(NULL);
  max_its = MAX_ITS;
  if (setup.max_iterations > 0) max_its = setup.max_iterations;

  /* to be safe, initialize overall potential to 0 */
  ig = lrint(gridstep[0]/grid);  // here ig is the number of grid-point distances for the coarsest grid we will use
  for (z=0; z<L+ig; z++) {
    for (r=0; r<R+ig; r++) {
      v[0][z][r] = v[1][z][r] = 0;
    }
  }
  /* ----------------------------------------------------------------------------
     now set up and perform the WP relaxation for each of the grid step sizes in turn
  */
  for (istep = 0; istep < 4 && gridstep[istep] > 0; istep++) {
    if (istep > 0) {
      /* not the first go-around, so the previous calculation was on a coarser grid...
	 now interpolate the potential to the final fine grid */
      interpolate(&setup, v, grid);
    }
    grid = gridstep[istep];
    old = 1;
    new = 0;
    ig = lrintf(grid / grid0);

    //-------------------------------------------------
    if (istep == 0) {
      // no previous coarse relaxation, so set initial potential:
      for (z=0; z<L+1; z++) {
	for (r=0; r<R+1; r++) {
	  v[0][z][r] = v[1][z][r] = 0.0;
	}
      }
      /*  ----- can comment out this next section to test convergence of WP ----- */
      // perhaps this is a better initial guess than just zero everywhere
      a = LC + RC / 2;
      b = 2.0f*a / (float) (L + R);
      for (z=1; z<L; z++) {
	for (r=1; r<R; r++) {
	  c = a / sqrt(z*z + r*r) - b;
	  if (c < 0) c = 0;
	  if (c > 1) c = 1;
	  v[0][z][r] = v[1][z][r] = c;
	}
      }
    }
    /* boundary conditions and permittivity
       boundary condition at Ge-vacuum interface:
       epsilon0 * E_vac = espilon_Ge * E_Ge
    */
    set_boundary_conditions(&setup, grid, v, bulk, vfraction, 0);
    if (setup.verbosity >= NORMAL) printf("grid = %f  ig = %d\n", grid, ig);
    //-------------------------------------------------
    /* In the following we divide areas and volumes by pi
       r_bin   rmax  A_top A_outside A_inside  volume  total_surf  out/top  tot/vol
       0     1/2    1/4      1         0       1/4      1.5         4        6  << special case
       1     3/2      2      3         1        2        8        3/2        4
       2     5/2      4      5         3        4       16        5/4        4
       3     7/2      6      7         5        6       24        7/6        4
       r   r+0.5     2r    2r+1      2r-1      2r       8r     (2r+1)/2r     4
       = 1+0.5/r
    */
    // weighting values for the relaxation alg. as a function of r
    s1[0] = 4.0;
    s2[0] = 0.0;
    for (r=ig; r<R+ig; r+=ig) {
      s1[r] = 1.0 + 0.5 / (double) (r/ig);   //  for r+1
      s2[r] = 1.0 - 0.5 / (double) (r/ig);   //  for r-1
    }
    for (z=0; z<L+ig; z+=ig) {
      for (r=0; r<R+ig; r+=ig) {
	eps[z][r] = eps_dz[z][r] = eps_dr[z][r] = 16;   // permittivity inside Ge
	if (z < LO  && r < RO && r > RO-WO-1) eps[z][r] =  1;  // permittivity inside vacuum
	if (r > 0) eps_dr[z][r-ig] = (eps[z][r-ig]+eps[z][r])/2.0f;
	if (z > 0) eps_dz[z-ig][r] = (eps[z-ig][r]+eps[z][r])/2.0f;
      }
    }

    // now do the actual relaxation
    for (z=0; z<L; z+=ig) {
      for (r=0; r<R; r+=ig) {
        d[z][r][0] = z-ig;
        d[z][r][1] = z+ig;
        d[z][r][2] = r-ig;
        d[z][r][3] = r+ig;
      }
    }
    if (fix_adaptive && istep > 1) {
      // see if the pixel is at least four grid steps away from a boundary (or r=0)
      // remember that ig is the curent pixel size (in units of grid0)
      for (z=4*ig; z<L-4*ig; z+=3*ig) {
        for (r=4*ig; r<R-4*ig; r+=3*ig) {
          if (bulk[z][r] || bulk[z][r-ig] ||
              bulk[z-ig][r] || bulk[z-ig][r-ig] || bulk[z-ig][r+ig]) continue;
          int ok = 1;
          for (i=-4*ig; i<=4*ig; i+=ig) {
            for (j=-4*ig; j<=4*ig; j+=ig) {
              if (bulk[z+i][r+j] &&
                  bulk[z+i][r+j]!= -4 &&
                  bulk[z+i][r+j]!=  4) ok=0; // no, it's not
            }
          }
          if (ok) {   // yes it is
            
            d[z][r][0] = z-2*ig;
            d[z][r][1] = z+2*ig;
            d[z][r][2] = r-2*ig;
            d[z][r][3] = r+2*ig;
            
            for (i=-ig; i<=ig; i++) {
              for (j=-ig; j<=ig; j++) {
                bulk[z+i][r+j] = -4;
              }
            }
            if (fix_adaptive > 2) bulk[z][r] = 4;
          }
        }
      }
    }
    int n=0;
    for (z=0; z<L; z+=ig) {
      for (r=0; r<R; r+=ig) {
        if (bulk[z][r] == -4 || bulk[z][r] == 4) n++;
      }
    }
    printf("\n ------ n = %d of %d\n", n, L*R/ig/ig);

    do_relax(R, L, LC, istep, max_its, v, bulk, ig, eps_dr, eps_dz,
             s1, s2, undepleted, imp, 0, d);

    if (setup.verbosity >= CHATTY) {
      t1 = time(NULL);
      printf(" ^^^^^^^^^^^^^ %d (%d) s elapsed ^^^^^^^^^^^^^^\n",
	     (int) (t1 - t0), (int) (t1 - t2));
      t2 = t1;
    }
    if (istep == 0) max_its /= MAX_ITS_FACTOR;
  }
  if (fix_adaptive == 2 || fix_adaptive == 4) {
    for (z=0; z<L; z++) {
      for (r=0; r<R; r++) {
        if (bulk[z][r] == -4 || bulk[z][r] == 4) {
          bulk[z][r] = 0;
          d[z][r][0] = z-ig;
          d[z][r][1] = z+ig;
          d[z][r][2] = r-ig;
          d[z][r][3] = r+ig;
        }
      }
    }   
    do_relax(R, L, LC, istep, 1000, v, bulk, ig, eps_dr, eps_dz,
             s1, s2, undepleted, imp, 0, d);
  }

  /* --------------------- calculate capacitance ---------------------
     1/2 * epsilon * integral(E^2) = 1/2 * C * V^2
     so    C = epsilon * integral(E^2) / V^2
     V = 1 volt
  */
  printf("Calculating integrals of weighting field\n");
  esum = esum2 = j = 0;
  for (z=0; z<L; z++) {
    for (r=1; r<R; r++) {
      E_r = eps_dr[z][r]/16.0 * (v[new][z][r] - v[new][z][r+1])/(0.1*grid);
      E_z = eps_dz[z][r]/16.0 * (v[new][z][r] - v[new][z+1][r])/(0.1*grid);
      esum += (E_r*E_r + E_z*E_z) * (double) r;
      if ((r == RC && z <= LC) ||
	  (r <= RC && z == LC) ||
	  (r == RC+1 && z <= LC+1) || // average over two different surfaces
	  (r <= RC+1 && z == LC+1)) {
	if (bulk[z+1][r+1] < 0) j = 1;
	esum2 += 0.5 * sqrt(E_r*E_r + E_z*E_z) * (double) r;  // 0.5 since averaging over 2 surfaces
      }
    }
  }
  esum  *= 2.0 * pi * 0.01 * Epsilon * pow(grid, 3.0);
  // Epsilon is in pF/mm
  // 0.01 converts (V/cm)^2 to (V/mm)^2, pow() converts to grid^3 to mm3
  esum2 *= 2.0 * pi * 0.1 * Epsilon * pow(grid, 2.0);
  // 0.1 converts (V/cm) to (V/mm),  grid^2 to  mm2
  printf("\n  >>  Calculated capacitance at %.0f V: %.3lf pF\n", BV, esum);
  if (j==0) {
    printf("  >>  Alternative calculation of capacitance: %.3lf pF\n\n", esum2);
  } else {
    printf("\n");
  }

  if (WP) {
    // write WP values to output file
    if (!(file = fopen(setup.wp_name, "w"))) {
      printf("ERROR: Cannot open file %s for weighting potential...\n", setup.wp_name);
      return 1;
    } else {
      printf("Writing weighting potential to file %s\n\n", setup.wp_name);
    }
    /* copy configuration parameters to output file */
    report_config(file, config_file_name);
    fprintf(file, "#\n# HV bias in fieldgen: %.1f V\n", BV);
    if (fully_depleted) {
      fprintf(file, "# Detector is fully depleted.\n");
    } else {
      fprintf(file, "# Detector is not fully depleted.\n");
      if (bubble_volts > 0.0f) fprintf(file, "# Pinch-off bubble at %.0f V potential\n", bubble_volts);
    }
    fprintf(file, "#\n## r (mm), z (mm), WP\n");
    if (WP > 1) {
      for (r=-R; r<R+1; r++) {
        if ((i=r) < 0) i = -r;
        for (z=0; z<L+1; z++) {
          fprintf(file, "%7.2f %7.2f %12.6e\n",
                  ((float) r)*grid,  ((float) z)*grid, v[new][z][i]);
        }
        fprintf(file, "\n");
      }
    } else {
      for (r=0; r<R+1; r++) {
        for (z=0; z<L+1; z++) {
          fprintf(file, "%7.2f %7.2f %12.6e\n",
                  ((float) r)*grid,  ((float) z)*grid, v[new][z][r]);
        }
        fprintf(file, "\n");
      }
    }
    fclose(file);
  }

  if (fully_depleted) {
    /* estimate depletion voltage */
    min = BV;
    for (z=0; z<LC+2; z++) {
      for (r=0; r<RC+2; r++) {
        if (vsave[z][r] > 0 &&
            min > vsave[z][r] / (1.0 - v[new][z][r])) {
          min = vsave[z][r] / (1.0 - v[new][z][r]);
          vminr = r;
          vminz = z;
        }
      }
    }

    /* also try to check for bubble depletion / pinch-off
       by seeing how much the bias must be reduced for any pixel to be in a local potential minimum */
    min2 = BV;
    vminr2 = vminz2 = 0;
    /* first check along r=0 (z-axiz) */
    for (z=LC+2; z < L - setup.hole_length/grid0 - 2; z++) {
      dVn[0] = vsave[z+1][0]  - vsave[z][0] ;
      dWn[0] = v[new][z+1][0] - v[new][z][0];
      dVn[1] = vsave[z-1][0]  - vsave[z][0] ;
      dWn[1] = v[new][z-1][0] - v[new][z][0];
      test = -1;
      for (i=0; i<2; i++) {
        if (dWn[i]*grid > 0.00001 && dVn[i] < 0 && test < -dVn[i]/dWn[i]) test = -dVn[i]/dWn[i];
      }
      if (test >= 0 && min2 > test) {
        min2 = test;
        vminr2 = 0;
        vminz2 = z;
      }
    }
    /* then check pinch-off for r > 0 */
    for (z=LC+2; z<L-2; z++) {
      for (r=1; r<R-2; r++) {
        if (vsave[z][r] > 0 && v[new][z][r] > 0.0001) {
          dVn[0] = vsave[z+1][r]  - vsave[z][r] ;
          dWn[0] = v[new][z+1][r] - v[new][z][r];
          dVn[1] = vsave[z-1][r]  - vsave[z][r] ;
          dWn[1] = v[new][z-1][r] - v[new][z][r];
          dVn[2] = vsave[z][r+1]  - vsave[z][r] ;
          dWn[2] = v[new][z][r+1] - v[new][z][r];
          dVn[3] = vsave[z][r-1]  - vsave[z][r] ;
          dWn[3] = v[new][z][r-1] - v[new][z][r];
          test = -1;
          for (i=0; i<4; i++) {
            if (dWn[i]*grid > 0.00001 && dVn[i] < 0 && test < -dVn[i]/dWn[i])
              test = -dVn[i]/dWn[i];
          }
          if (test >= 0 && min2 > test) {
            min2 = test;
            vminr2 = r;
            vminz2 = z;
          }
        }
      }
    }
    if (min2 < min) {
      printf("\nEstimated pinch-off voltage = %.0f V\n", BV - min);
      printf(" min2 = %.1f at (r,z) = (%.1f, %.1f), so\n", min2, vminr2*grid, vminz2*grid);
      printf("   Full depletion (max pinch-off voltage) = %.0f\n", BV - min2);
    } else {
      printf("\nEstimated depletion voltage = %.0f V\n", BV - min);
    }
    printf("\nMinimum bulk field = %.1f V/cm at (r,z) = (%.1f, %.1f) mm\n\n",
           Emin, rmin, zmin);
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

#define SQ(x) ((x)*(x))
int outside_detector(float z, float r, MJD_Siggen_Setup *setup) {
  float lith  = setup->Li_thickness;
  float zmax  = setup->xtal_length - lith;
  float rmax  = setup->xtal_radius - lith;
  float r1, z1, br, a, b;

  if (z >= zmax || z <= 0) return 1;
  if (r < 0) r = -r;
  if (r >= rmax) return 1;
  if (r >= setup->wrap_around_radius && z <= lith) return 1;

  r1 = rmax - r;  // distance from outer radius
  z1 = zmax - z;  // distance from top of crystal
  /* check point contact */
  if (z <= setup->pc_length && r <= setup->pc_radius) {
    if (!setup->bulletize_PC) return 1;
    if (setup->pc_length > setup->pc_radius) {
      a = setup->pc_length - setup->pc_radius;
      if (z <= a || SQ(z-a) + SQ(r) <= SQ(setup->pc_radius)) return 1;
    } else {
      a = setup->pc_radius - setup->pc_length;
      if (r <= a || SQ(z) + SQ(r-a) <= SQ(setup->pc_length)) return 1;
    }
    return 0;
  }
  /* check ditch */
  if (z <= setup->ditch_depth  &&
      setup->ditch_thickness > 0 && setup->wrap_around_radius > 0 &&
      r <= setup->wrap_around_radius &&
      r >= setup->wrap_around_radius - setup->ditch_thickness) return 1;

  /* check hole */
  if ( r <= setup->hole_radius + lith &&
      z1 <= setup->hole_length) {  // note no lith added here since it was already subtracted from zmax
    br = setup->hole_bullet_radius + 0.8*lith;
    b = zmax - setup->hole_length + br;
    if (z >= b) return 1;
    a = setup->hole_radius + lith - br;
    if (r <= a || SQ(b-z) + SQ(a-r) <= SQ(br)) return 1;
  }
  /* check inner taper of hole */
  if (z1 <= setup->inner_taper_length &&
      r  <= setup->hole_radius + lith +
            ((setup->inner_taper_length - z1) *
              setup->inner_taper_width / setup->inner_taper_length)) return 1;      
  /* check outer taper of crystal */
  if (z1 <= setup->outer_taper_length &&
      r1 <= ((setup->outer_taper_length - z1) *
              setup->outer_taper_width / setup->outer_taper_length)) return 1;
  /* check 45-degree bottom outer taper of crystal */
  if ( z <= setup->bottom_taper_length + 0.71*lith &&
      r1 <= setup->bottom_taper_length + 0.71*lith - z) return 1;

  /* check bulletizations */
  br = setup->top_bullet_radius - lith * 0.7;
  if (z1 <= br &&
      r1 <= br &&
      SQ(br - r1) + SQ(br - z1) >= br*br) return 1;
  br = setup->bottom_bullet_radius - lith * 0.7;
  if ( z <= br + lith &&
      r1 <= br &&
      SQ(br - r1) + SQ(br - z + lith) >= br*br) return 1;

  return 0;
}
#undef SQ

int set_boundary_conditions(MJD_Siggen_Setup *setup, float grid, double **v[2],
                            int **bulk, double **vfraction, float BV)
{

  /* NOTE: BV is bias voltage for EV calculation, or 0 for WP calculation */
  float  grid0 = setup->xtal_grid;  // final desired grid size
  int    L = lrint(setup->xtal_length/grid0);
  int    R = lrint(setup->xtal_radius/grid0);
  int    ig, z, r, n;
  float  fz, fr, d = grid /2.1;

  ig = lrintf(grid / grid0);

  for (z=0; z<L+ig; z++) {
    for (r=0; r<R+ig; r++) {
      bulk[z][r] = 0;
      vfraction[z][r] = 1.0;
    }
  }
  for (z=0; z<L+ig; z+=ig) {
    for (r=0; r<R+ig; r+=ig) {
      fz = grid0 * z;
      fr = grid0 * r;
      n = 0;
      n += outside_detector(fz+d, fr-d, setup);
      n += outside_detector(fz+d, fr+d, setup);
      n += outside_detector(fz-d, fr-d, setup);
      n += outside_detector(fz-d, fr+d, setup);

      if (n == 0) continue; // no complications

      n += outside_detector(fz, fr, setup);
      n += outside_detector(fz, fr-d, setup);
      n += outside_detector(fz, fr+d, setup);
      n += outside_detector(fz-d, fr, setup);
      n += outside_detector(fz+d, fr, setup);
      vfraction[z][r] = 1.0 - n/9.0;
      if (n > 4) {
        if (fz-d > setup->ditch_depth ||
            fr > setup->wrap_around_radius) {
          bulk[z][r] = -1;           // fixed contact voltage
          v[0][z][r] = v[1][z][r] = BV;
        } else if (fz-d < setup->pc_length &&
                   fr-d < setup->pc_radius) {
          bulk[z][r] = -1;           // fixed contact voltage
          v[0][z][r] = v[1][z][r] = (BV != 0? 0 : 1);
        } else {
          bulk[z][r] = 3;            // ditch
        }
      }
    }
  }

  /* for pixels adjacent to the ditch, we also need to set bulk to 3 */
  for (z=0; z<L; z+=ig) {
    for (r=0; r<R; r+=ig) {
      if (bulk[z][r] != 3) continue;
      if (z > 0 && bulk[z-ig][r] == 0) bulk[z-ig][r] = -3;
      if (bulk[z+ig][r] == 0) bulk[z+ig][r] = -3;
      if (r > 0 && bulk[z][r-ig] == 0) bulk[z][r-ig] = -3;
      if (bulk[z][r+ig] == 0) bulk[z][r+ig] = -3;
    }
  }
  for (z=0; z<L+ig; z+=ig) {
    for (r=0; r<R+ig; r+=ig) {
      if (bulk[z][r] == -3) bulk[z][r] = 3;
    }
    if (bulk[z][0] == 0) bulk[z][0] = 3;
  }
  /* for pixels at z=0 or r=0, we also need to set bulk to 3 */
  for (z=0; z<L+ig; z+=ig)
    if (bulk[z][0] == 0) bulk[z][0] = 3;
  for (r=0; r<R+ig; r+=ig)
    if (bulk[0][r] == 0) bulk[0][r] = 3;

  return 0;   
}

int do_relax(int R, int L, int LC, int istep, int max_its, double **v[2], int **bulk, int ig,
             double **eps_dr, double **eps_dz, double *s1, double *s2,
             char **undepleted, float **imp, int ev_calc, int **d[6]) {

  int    old = 1, new = 0, iter, r, z, *p;
  double eps_sum, v_sum, mean, save_dif;
  float  dif, sum_dif, max_dif;
  int    adaptive = 0;


// the following definition of the factor for over-relaxation improves convergence time by a factor ~ 70 for a 2kg ICPC detector
#define OVER_RELAX_FACTOR (0.988 * (1.0 - 1.0/(double)(1+iter/4)))
  //#define OVER_RELAX_FACTOR   (0.989 - 0.0015*ig)   // CHECKME

  for (iter=0; iter<max_its; iter++) {
    double OR_fact = OVER_RELAX_FACTOR;
    if (iter < 2) OR_fact = 0.0;
    else if (iter < 200) OR_fact *= 0.95;   // CHECKME

    old = 1 - old;
    new = 1 - new;
    sum_dif = 0.0f;
    max_dif = 0.0f;

    for (z=0; z<L; z+=ig) {
      for (r=0; r<R; r+=ig) {
        if (bulk[z][r] < 0) continue;      // outside or inside contact, or fixed adaptive grid point
        save_dif = v[old][z][r] - v[new][z][r];  // step difference from previous iteration
        p = d[z][r];
        /* p[] gives z or r position of neighbor pixels, e.g.:
           p[0] = d[z][r][0] = z-ig;
           p[1] = d[z][r][1] = z+ig;
           p[2] = d[z][r][2] = r-ig;
           p[3] = d[z][r][3] = r+ig;
        */

        if (bulk[z][r] == 0) {            // normal bulk, no complications
          v_sum = (v[old][p[1]][r] + v[old][z][p[3]]*s1[r] +
                   v[old][p[0]][r] + v[old][z][p[2]]*s2[r]);
          eps_sum = 4;
        } else if (bulk[z][r] == 4) {            // larger pixel (3x larger on each side)
          v_sum = (v[old][p[1]][r] + v[old][z][p[3]]*(1.0 + 1.5 / (double) (r / ig)) +
                   v[old][p[0]][r] + v[old][z][p[2]]*(1.0 - 1.5 / (double) (r / ig)));
          eps_sum = 4;
        } else if (bulk[z][r] == 3) {            // in or adjacent to the ditch, or at z = 0 or r=0
          v_sum = v[old][p[1]][r]*eps_dz[z][r] + v[old][z][p[3]]*eps_dr[z][r]*s1[r];
          eps_sum = eps_dz[z][r] + eps_dr[z][r]*s1[r];
          if (z > 0) {
            v_sum += v[old][p[0]][r]*eps_dz[p[0]][r];
            eps_sum += eps_dz[p[0]][r];
          } else {
            v_sum += v[old][p[1]][r]*eps_dz[z][r];  // reflection symm around z=0
            eps_sum += eps_dz[z][r];
          }
          if (r > 0) {
            v_sum += v[old][z][p[2]]*eps_dr[z][p[2]]*s2[r];
            eps_sum += eps_dr[z][p[2]]*s2[r];
          } else {
            v_sum += v[old][z][p[3]]*eps_dr[z][r]*s1[r];  // reflection symm around r=0
            eps_sum += eps_dr[z][r]*s1[r];
          }
        } else if (bulk[z][r] == 1) {    // interpolated radial edge of point contact
          /* since the PC radius is not in the middle of a pixel,
             use a modified weight for the interpolation to (p[2])
          */
          v_sum = v[old][p[1]][r]*eps_dz[z][r] + v[old][z][p[3]]*eps_dr[z][r]*s1[r] +
            v[old][z][p[2]]*eps_dr[z][p[2]]*s2[r];//*frrc[z];
          eps_sum = eps_dz[z][r] + eps_dr[z][r]*s1[r] + eps_dr[z][p[2]]*s2[r];//*frrc[z];
          if (z > 0) {
            v_sum += v[old][p[0]][r]*eps_dz[p[0]][r];
            eps_sum += eps_dz[p[0]][r];
          } else {
            v_sum += v[old][p[1]][r]*eps_dz[z][r];  // reflection symm around z=0
            eps_sum += eps_dz[z][r];
          }
        } else if (bulk[z][r] == 2) {    // interpolated z edge of point contact
          /* since the PC length is not in the middle of a pixel,
             use a modified weight for the interpolation to (p[0])
          */
          v_sum = v[old][p[1]][r]*eps_dz[z][r] + v[old][z][p[3]]*eps_dr[z][r]*s1[r] +
            v[old][p[0]][r]*eps_dz[p[0]][r];//*fLC;
          eps_sum = eps_dz[z][r] + eps_dr[z][r]*s1[r] + eps_dz[p[0]][r];//*fLC;
          if (r > 0) {
            v_sum += v[old][z][p[2]]*eps_dr[z][p[2]]*s2[r];
            eps_sum += eps_dr[z][p[2]]*s2[r];
          } else {
            v_sum += v[old][z][p[3]]*eps_dr[z][r]*s1[r];  // reflection symm around r=0
            eps_sum += eps_dr[z][r]*s1[r];
          }
          // check for cases where the PC corner needs modification in both r and z
          if (z == LC && bulk[p[0]][r] == 1) {
            v_sum += v[old][z][p[2]]*eps_dr[z][p[2]]*s2[r];//*(frrc[z]-1.0);
            eps_sum += eps_dr[z][p[2]]*s2[r];//*(frrc[z]-1.0);
          }

        } else {
          printf(" ERROR! bulk = %d undefined for (z,r) = (%d,%d)\n",
                 bulk[z][r], z, r);
          return 1;
        }

        // calculate the interpolated mean potential and the effect of the space charge
        mean = v_sum / eps_sum;
        v[new][z][r] = mean + ev_calc * imp[z][r];

        // calculate difference from last iteration, for convergence check
        dif = v[old][z][r] - v[new][z][r];
        if (dif < 0) dif = -dif;
        sum_dif += dif;
        if (max_dif < dif) max_dif = dif;
        // do over-relaxation
        if (adaptive)
          //v[new][z][r] += 0.898 * OR_fact*save_dif;
          v[new][z][r] += (0.966 - 0.07*ev_calc) * OR_fact*save_dif;
        else
          v[new][z][r] += OR_fact*save_dif;

        //if (v[new][z][r] < 0) v[new][z][r] = 0;

        // fill in changes for neighboring pixels for large adaptive-grid pixels
        if (bulk[z][r] == 4) {
          adaptive = 1;
          dif = v[new][z][r] - v[old][z][r];
          // dif *= 1.01;
          int i, j;
          for (i=-ig; i<=ig; i++) {
            for (j=-ig; j<=ig; j++) {
              v[new][z+i][r+j] = v[old][z+i][r+j] + dif;
            }
          }
          if (0 && iter == 1900) {
            printf("z,r = %d %d; dif = %.2e, v = ", z, r, dif);
            for (i=-ig; i<=ig; i++) {
              for (j=-ig; j<=ig; j++) {
                printf(" %.2e", v[new][z+i][r+j]);
              }
            }
            printf("\n");
          }
        }

      }
    }

    // report results for some iterations
    if (iter < 10 || (iter < 600 && iter%100 == 0) || iter%1000 == 0) {
      if (ev_calc) {
        printf("%5d %d %d %.10f %.10f\n", iter, old, new, max_dif, sum_dif/(float) (L*R));
      } else {
        printf("%5d %d %d %.10f %.10f ; %.10f %.10f\n",
               iter, old, new, max_dif, sum_dif/(float) (L*R),
               v[new][L/2][R/2], v[new][L-5][R-5]);
      }
    }
    // check for convergence
    //if ( ev_calc && max_dif < 0.00000008/ig) break;
    if ( ev_calc && max_dif < 0.00000008) break;
    if (!ev_calc && max_dif < 0.0000000001) break;

    if (adaptive && iter%50 == 0)
      if (ev_calc) interpolate2(v, bulk, ig*3, L, R);

  }
  printf(">> %d %.16f\n\n", iter, sum_dif);
  if (adaptive && ev_calc) interpolate2(v, bulk, ig, L, R);

  return 0;
}

int interpolate(MJD_Siggen_Setup *setup, double **v[2], float grid) {

  double f, f1z, f2z, f1r, f2r;
  float  grid0 = setup->xtal_grid;  // final desired grid size
  int    L = lrint(setup->xtal_length/grid0);
  int    R = lrint(setup->xtal_radius/grid0);
  int    ig, z, r, zz, rr, rmax, zmax;

  ig = lrintf(grid / grid0);
  if (ig < 2) {
    printf("Interpolate(): ig < 2; grid = %f  grid0 =%f\n", grid, grid0);
    return 0;
  }
  f = 1.0 / (float) ig;
  printf("\n Interpolating grid %.4f -> %.4f  f = %f\n\n", grid, grid0, f);
  for (z=0; z<L+1; z+=ig) {
    for (r=0; r<R+1; r+=ig) {
      f1z = 0.0;
      zmax = z+ig;
      if (zmax > L+1) zmax = L+1;
      for (zz=z; zz<z+ig; zz++) {
        f2z = 1.0 - f1z;
        f1r = 0.0;
        rmax = r+ig;
        if (rmax > R+1) rmax = R+1;
        for (rr=r; rr<rmax; rr++) {
          f2r = 1.0 - f1r;
          v[0][zz][rr] =      // linear interpolation of potential
            f2z*f2r*v[1][z][r   ] + f1z*f2r*v[1][z+ig][r   ] +
            f2z*f1r*v[1][z][r+ig] + f1z*f1r*v[1][z+ig][r+ig];
          /*
            DCR Nov 2019
            If the potential here is sufficiently linear,
            and the grid size is a small fraction of the crystal size,
            and the detector is fully depleted,
            then we can fix the potential to the interpolated value from the previous istep
            so set bulk to special value of -2
          */
          /*
          if (fix_adaptive && fully_depleted &&
              gridstep[istep] < setup.xtal_radius/80.0 &&
              gridstep[istep] < setup.xtal_length/80.0 &&
              z > 1 && z < L-1 && r > 1 && r < R-1 && v[1][z][r] !=0 &&
              fabs((v[1][z][r-1] + v[1][z][r+1] - 2.0*v[1][z][r]) / v[1][z][r]) < 0.0005 &&
              fabs((v[1][z-1][r] + v[1][z+1][r] - 2.0*v[1][z][r]) / v[1][z][r]) < 0.0005) {
          }
          */
          f1r += f;
        }
        f1z += f;
      }
    }
  }
  /* be sure to copy new potential map to the other (old=1) version */
  for (z=0; z<L+1; z++) {
    for (r=0; r<R+1; r++) {
      v[1][z][r] = v[0][z][r];
    }
  }

  return 0;
}

int interpolate2(double **v[2], int **bulk, int ig, int L, int R) {

  double f, ez, er, ezr1, erz1, ezr2, erz2;
  int    ng, z, r, zz, rr, z1, r1, z2, r2;

  // printf("Interpolate2:  ig = %d\n", ig);
  if (ig < 2) return 0;
  f = 0.5 / (float) ig;
  ng = ig/2;
  for (z=0; z<L+1; z++) {
    z1 = z - ig;
    z2 = z + ig;
    for (r=0; r<R+1; r++) {
      if (bulk[z][r] != 4) continue;
      r1 = r - ig;
      r2 = r + ig;
      ez   = f   * (v[1][z1][r ] - v[1][z2][r ]);
      ezr1 = f*f * (v[1][z1][r1] - v[1][z2][r1]) + (1.0-f) * ez;
      ezr2 = f*f * (v[1][z1][r2] - v[1][z2][r2]) + (1.0-f) * ez;
      er   = f   * (v[1][z ][r1] - v[1][z ][r2]);
      erz1 = f*f * (v[1][z1][r1] - v[1][z1][r2]) + (1.0-f) * er;
      erz2 = f*f * (v[1][z2][r1] - v[1][z2][r2]) + (1.0-f) * er;
      for (zz=z-ng; zz<=z+ng; zz++) {
        for (rr=r-ng; rr<=r+ng; rr++) {
          if (zz < z) {
            if (rr < r) {
              v[0][zz][rr] = v[1][zz][rr] = v[1][z][r] + ezr1 + erz1;
            } else if (rr > r) {
              v[0][zz][rr] = v[1][zz][rr] = v[1][z][r] + ezr2 - erz1;
            } else {
              v[0][zz][rr] = v[1][zz][rr] = v[1][z][r] + ez;
            }
          } else if (zz > z) {
            if (rr < r) {
              v[0][zz][rr] = v[1][zz][rr] = v[1][z][r] - ezr1 + erz2;
            } else if (rr > r) {
              v[0][zz][rr] = v[1][zz][rr] = v[1][z][r] - ezr2 - erz2;
            } else {
              v[0][zz][rr] = v[1][zz][rr] = v[1][z][r] - ez;
            }
          } else {
            if (rr < r) {
              v[0][zz][rr] = v[1][zz][rr] = v[1][z][r] + er;
            } else if (rr > r) {
              v[0][zz][rr] = v[1][zz][rr] = v[1][z][r] - er;
            }
          }
        }
      }
    }
  }

  return 0;
}
