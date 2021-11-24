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
  int    i, j, t;
  float  e, r, z, *s;
  int    his[3][4096] = {{0}};
  point  pt;
  double n0=0, n1=0, n2=0, n3=0, s1, s2, s3, gvol;
  FILE  *fp;

  
  /* read config file and initialize */
  if (argc < 2) {
    printf("Usage: %s <config_file_name>\n", argv[0]);
    return 1;
  }
  if (signal_calc_init(argv[1], &setup) != 0) return 1;

  setup.preamp_tau       = 30;
  setup.step_time_calc   = 5;
  setup.step_time_out    = 10;
  setup.time_steps_calc  = 6000;
  setup.ntsteps_out      = 3000;
  setup.verbosity        = 0;
  gvol = 3.1416/pow(setup.xtal_grid, 3.0);
  if ((s = malloc(setup.ntsteps_out*sizeof(*s))) == NULL) {
    printf("malloc failed\n");
    return 1;
  }
  
  for (r = 0; r < setup.rlen; r += 1.0) {
    // printf("r = %.1f\n", r);
    for (z = 0; z < setup.zlen; z += 1.0) {
      pt.x = r;
      pt.y = 0;
      pt.z = z;

      if (get_signal(pt, s, &setup) > 0) {
        s1 = s2 = s3 = 0;
        for (t = 0; t < setup.ntsteps_out; t++){
          if (s[t] > 0.99) break;
          s3 += 1.0-s[t];
          if (s[t] < 0.001) s1 = s3;
          if (s[t] < 0.004) s2 = s3;
        }
        his[0][lrintf(s1)]++;
        his[1][lrintf(s2)]++;
        his[2][lrintf(s3)]++;
      }

      n0 += r*r*gvol;
      e = setup.wpot[lrint(r)][lrintf(z)];
      if (e > 1e-12) n1 += r*r*gvol;
      if (e > 1e-3)  n2 += r*r*gvol;
      if (e > 4e-3)  n3 += r*r*gvol;
    }
  }

  // n0 = setup.rlen * (setup.rlen + 1) / 2  * setup.zlen;
  printf("\n  ---- WP: ----\n");
  printf ("Total grid-point volume : %6.1f\n", n0);
  printf ("        With nonzero WP : %6.1f (%4.1f%% of total)\n",        n1, 100.0 * n1 / n0);
  printf ("         Threshold 0.1%% : %6.1f (%4.1f%% of non-zero WP)\n", n2, 100.0 * n2 / n1);
  printf ("         Threshold 0.4%% : %6.1f (%4.1f%% of non-zero WP)\n", n3, 100.0 * n3 / n1);

  for (i = 0; i < 3; i++) {
    s1 = s2 = s3 = 0;
    for (j = 0; j < 800; j++) {
      s1 += his[i][j];
      s2 += his[i][j] * j;
      s3 += his[i][j] * j*j;
    }
    s2 /= s1;
    s3 /= s1;
    t = s2;  // centroid;
    s3 -= s2*s2;
    s3 = 2.355*sqrt(s3)/20.0; // fwhm / 20
    for (j = 0; j < 800; j++) his[i][2000+t-j] = his[i][j];
    printf(" FWHM%d = %.2f (%.2f; %.2f)\n", i, s3,
           sqrt(s3*s3 + 2.25*2.25),  sqrt(s3*s3/4.0 + 2.25*2.25));
  }

  fp = fopen("dtc.spn", "w");
  fwrite(his[0], sizeof(his[0]), 3, fp);
  fclose(fp);

  return 0;
}
