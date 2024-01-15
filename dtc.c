#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "mjd_siggen.h"
#include "calc_signal.h"
#include "detector_geometry.h"
#include "fields.h"


int main(int argc, char **argv) {

  MJD_Siggen_Setup setup;
  int    i, j, t, tmax = 0, do_diffusion = 0;
  float  e, r, z, *s, value;
  int    his[4][4096] = {{0}};
  point  pt;
  double n0=0, n1=0, n2=0, n3=0, s1, s2, s3, gvol, a, amax;
  char   *description[3] = {"0.1% threshold", "0.4% threshold", "No"};
  FILE  *fp;
  // float  factor = 0.04, fwhm0 = 2.3;      // initial guesses
  float  factor = 0.0188, fwhm0 = 2.3325;       // calibrated from Mirion HADES data
  // float  factor = 0.0222, fwhm0 = 2.47;  // calibrated from ORTEC acceptance data
  // float  factor = 0.083, fwhm0 = 2.13;   // calibrated from ORTEC MJD data

  
  /* read config file and initialize */
  if (argc < 2) {
    printf("Usage: %s <config_file_name>\n", argv[0]);
    return 1;
  }
  if (strstr(argv[1], "-dd")) {
      do_diffusion = 1;
      argv[1] = argv[2];
  }
  if (signal_calc_init(argv[1], &setup) != 0) return 1;

  // handle same command line pararmeters at fieldgen geometry optimization code...
  for (i=2+do_diffusion; i<argc; i++) {
    if (strstr(argv[i], "-b") ||
        strstr(argv[i], "-w") ||
        strstr(argv[i], "-p") ||
        strstr(argv[i], "-r") ||
        strstr(argv[i], "-z") ||
        strstr(argv[i], "-o")) {
      i++; continue;

    } else if (strstr(argv[i], "-dd")) {
      do_diffusion = 1;
      printf(">>> Difussion , charge cloud size, self-repulsion added\n"); 
    } else if (strstr(argv[i], "-d")) {
      continue;
    } else if (strstr(argv[i], "-g")) {
      printf("argv[%d] = %s\n", i, argv[i]);
      if (*(argv[i]+2) >= 'a' && *(argv[i]+2) <= 'z') {
        j = -1;
        if (*(argv[i]+2) == 'w') j = 0; //   wrap_around_radius
        if (*(argv[i]+2) == 'g') j = 1; //   hole_length_gap
        if (*(argv[i]+2) == 'h') j = 2; //   hole_radius
        if (*(argv[i]+2) == 't') j = 3; //   inner_taper_length
        if (*(argv[i]+2) == 'a') j = 4; //   taper_angle
        if (*(argv[i]+2) == 'l') j = 5; //   xtal length
        if (*(argv[i]+2) == 'r') j = 6; //   xtal radius
        if (*(argv[i]+2) == 'z') j = 7; //   z-cut position in mm (adds to any other -z input)
      } else {
        j = atoi(argv[i]+2);
      }
      value = atof(argv[++i]);
      if (j==0) {
        setup.wrap_around_radius = value;
        printf(" g%d override: wrap_around_radius = %.1f\n", j, value);
      } else if (j==1) {
        setup.hole_length = setup.xtal_length - value;
        printf(" g%d override: hole_length_gap = %.1f, hole_length = %.1f\n", j, value, setup.hole_length);
        if (setup.inner_taper_length > setup.hole_length) {
          setup.inner_taper_length = setup.hole_length;
          printf("              ... and inner_taper_length = %.1f\n", setup.inner_taper_length);
        }
      } else if (j==2) {
        setup.hole_radius = value;
        printf(" g%d override: hole_radius = %.1f\n", j, value);
      } else if (j==3) {
        setup.inner_taper_length = value;
        printf(" g%d override: inner_taper_length = %.1f\n", j, value);
      } else if (j==4) {
        setup.taper_angle = value;
        printf(" g%d override: taper_angle = %.1f\n", j, value);
      } else if (j==5) {
        printf(" g%d override: xtal_length = %.1f\n", j, value);
        if (setup.xtal_length != value) {
          // maintain hole_length_gap by adjusting hole_length
          setup.hole_length -= setup.xtal_length - value;
          printf("              ... and hole_length = %.1f\n", setup.hole_length);
          if (setup.inner_taper_length > setup.hole_length) {
            setup.inner_taper_length = setup.hole_length;
            printf("              ... and inner_taper_length = %.1f\n", setup.inner_taper_length);
          }
        }
        setup.xtal_length = value;
      } else if (j==6) {
        setup.xtal_radius = value;
        printf(" g%d override: xtal_radius = %.1f\n", j, value);
      } else if (j==7) {
      } else {
        printf("\nERROR: illegal geometry parameter override %s\n\n", argv[i-1]);
        return -1;
      }
    } else {
      printf("Possible options:\n"
             "      -g<n> <value>   override geometry spec from config file with <value>\n"
             "          n = 0 : wrap_around_radius;   n = 1: hole_length_gap;     n = 2: hole_radius\n"
             "          n = 3 : inner_taper_length;   n = 4: taper_angle\n"
             "          n = 5 : xtal_length;          n = 6: xtal_radius\n"
             "      ignored: -b -w -p -d -r -z\n");
      return 1;
    }
  }
  if (setup.inner_taper_length > setup.hole_length) {
    setup.inner_taper_length = setup.hole_length;
    printf(" Warning: inner_taper_length limited to hole_length = %.1f\n", setup.hole_length);
  }
  // set some siggen parameters to standard values
  if (field_setup(&setup) != 0) return -1;

  setup.preamp_tau       = 30;
  setup.step_time_calc   = 5;
  setup.step_time_out    = 10;
  setup.time_steps_calc  = 6000;
  setup.ntsteps_out      = 3000;
  setup.verbosity        = 0;

  if (do_diffusion) {
    setup.charge_cloud_size = 1.5; setup.use_diffusion = 1; setup.energy = 2000;
  } else {
    setup.charge_cloud_size = 0;   setup.use_diffusion = 0; setup.energy = 0;
  }
  
  gvol = 2.0*3.1416/pow(setup.xtal_grid, 3.0);
  if ((s = malloc(setup.ntsteps_out*sizeof(*s))) == NULL) {
    printf("malloc failed\n");
    return 1;
  }
  printf("size, difusion, energy: %f %d %f\n",
         setup.charge_cloud_size, setup.use_diffusion, setup.energy);

  // loop over the detector volume in steps of 1 mm
  // for (r = 0; r < setup.rlen; r += 1.0) {
  for (r = 0; r < setup.rmax; r += 1.0) {
    // printf("r = %.1f\n", r);
    // for (z = 0; z < setup.zlen; z += 1.0) {
    for (z = 0; z < setup.zmax; z += 1.0) {
      pt.x = r;
      pt.y = 0;
      pt.z = z;

      if (get_signal(pt, s, &setup) > 0) {
        // if (r == 35.0) printf("z = %.0f\n", z);
        s1 = s2 = s3 = 0;
        for (t = 0; t < setup.ntsteps_out; t++){
          if (s[t] > 0.99) break;
          s3 += 1.0-s[t];
          if (s[t] < 0.001) s1 = s3;
          if (s[t] < 0.004) s2 = s3;
        }
        if (t > setup.ntsteps_out-1) {
          // printf("no sig at (r,z) = (%.0f,%.0f)\n", r, z);
          continue;
        }
        if (tmax < t) tmax = t;

        // calculate A/E
        amax = a = (s[3]+s[4]+s[5]) - (s[0]+s[1]+s[2]);
        for (i=4; i<t; i++) {
          a += s[i+2] - s[i-1] - s[i-2] + s[i-4];
          if (amax < a) amax = a;
        }

        // increment by r to really sample volume
        his[0][lrintf(s1)] += (int) (r+0.5);    // Qdrift distribution for signal threshold < 0.1%
        his[1][lrintf(s2)] += (int) (r+0.5);    // Qdrift distribution for signal threshold < 0.4%
        his[2][lrintf(s3)] += (int) (r+0.5);    // full Qdrift distribution
        his[3][lrintf(amax*400.0)] += (int) (r+0.5);    // A/E distribution
      }

      n0 += r*gvol;
      e = setup.wpot[lrint(r/setup.xtal_grid)][lrintf(z/setup.xtal_grid)];
      if (e > 1.0e-12) n1 += r*gvol;
      if (e > 1.0e-3)  n2 += r*gvol;
      if (e > 4.0e-3)  n3 += r*gvol;
    }
  }
  tmax *= (int) setup.step_time_out;

  // n0 = setup.rlen * (setup.rlen + 1) / 2  * setup.zlen;
  printf("\n  ---- WP: ----\n");
  printf ("Total grid-point volume : %12.0f\n", n0);
  printf ("        With nonzero WP : %12.0f (%4.1f%% of total)\n",        n1, 100.0 * n1 / n0);
  printf ("         0.1%% threshold : %12.0f (%4.1f%% of non-zero WP)\n", n2, 100.0 * n2 / n1);
  printf ("         0.4%% threshold : %12.0f (%4.1f%% of non-zero WP)\n", n3, 100.0 * n3 / n1);

  printf(" Maximum drift time = %d ns\n"
         " FWHM values (keV) at 2615 keV (assuming %.3f keV with no trapping):\n"
         "   added    bad trapping   normal trapping   assumption\n", tmax, fwhm0);
  fp = fopen("dtc.spn", "w");
  fwrite(his[0], sizeof(his[0]), 3, fp);
  fclose(fp);
  for (i = 0; i < 3; i++) {
    s1 = s2 = s3 = 0;
    for (j = 0; j < 800; j++) {
      s1 += his[i][j];
      s2 += his[i][j] * j;
      s3 += his[i][j] * j*j;
    }
    printf("s1, s2, s3: %.0lf %.0lf %.0lf\n", s1, s2, s3); fflush(stdout);
    s2 /= s1;
    s3 /= s1;
    t = s2;  // centroid;
    s3 -= s2*s2;  // variance
    s3 = 2.355*sqrt(s3) * factor;
    for (j = 0; j < 800; j++) his[i][2000+t-j] = his[i][j];
    printf("  %6.3f %10.3f %15.3f        %s DTC\n", s3, sqrt(s3*s3*3.0 + fwhm0*fwhm0),
             sqrt(s3*s3 + fwhm0*fwhm0), description[i]);
  }
  s1 = s2 = s3 = 0;
  for (i = 0; i < 1000; i++) {
    s1 += his[3][i];
    if (s2 < his[3][i]) {
      s2 = his[3][i];
      j = i;
    }
  }
  s2 = s3 = 0;
  for (i = j-10; i<j+10; i++) {
    s2 += his[3][i];
    s3 += his[3][i] * (float) i;
  }
  j = (int) (s3/s2 + 0.5);   // centroid of main peak

  s2 = s3 = 0;
  for (i = 0;   i < j-10; i++) s2 += his[3][i];
  for (i = j+8; i < 1000; i++) s3 += his[3][i];
  printf("j, s2, s3: %d %.0lf %.0lf\n", j, s2, s3); fflush(stdout);
  s2 *= 100.0/s1;
  s3 *= 100.0/s1;
  printf(" %6.3f%%  %6.3f%%  %6.3f%%  A/E fraction\n", s2, s3, s2+s3);

  fp = fopen("dtc.spn", "w");
  fwrite(his[0], sizeof(his[0]), 4, fp);
  fclose(fp);

  return 0;
}
