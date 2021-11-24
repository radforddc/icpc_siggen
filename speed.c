/*  
DCR: quick and ugly hack to get total speed of holes as a function of position
     for plotting purposes. Just used for PPC/BEGe detectors.

 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>

#include "mjd_siggen.h"
#include "calc_signal.h"
#include "cyl_point.h"
#include "detector_geometry.h"
#include "fields.h"


int main(int argc, char **argv) {

  MJD_Siggen_Setup setup;

  int   i, j, n, seg, vel=1;
  int   ni[25] = {0};
  float *s;
  float ri[25][1000], zi[25][1000], lneg, lpos;
  float new, old = 25000.0, f;
  float r, phi, z, r45 = 3.1416/4.0;

  struct cyl_pt cyl;
  struct point cart;
  point *dpath_h;
  vector v;
  FILE  *file;

  int pad = 1;  // set to 1 to write zero values for plotting; required for python

  
  /* read config file and initialize */
  if (argc < 2) {
    printf("Usage: %s <config_file_name> [-t]\n"
           "  Specify -t to calculate drift times and isochrones instead of drift speed.\n",
           argv[0]);
    return 1;
  }
  if (argc > 2 && strstr(argv[2], "-t")) vel = 0;

  if (signal_calc_init(argv[1], &setup) != 0) return 1;
  setup.coord_type = CYL;

  if (vel) {
    file = fopen("speed.dat", "w");
    fprintf(file, "## r (mm), z (mm), Speed (cm/microsec)\n");
  } else {
    file = fopen("time.dat", "w");
    fprintf(file, "## r (mm), z (mm), drift-time (ns)\n");
  }

  if ((s = malloc(setup.time_steps_calc*sizeof(*s))) == NULL) {
    printf("malloc failed\n");
    return 1;
  }

  for (r = -setup.xtal_radius; r <= setup.xtal_radius; r+=0.25) {
    printf("r = %.2f\n", r);
    old = -1; // set negative to indicate that we should not yet look for isochrones
    for (z = 0.25; z <= setup.xtal_length; z+=0.25) {
      if (r < -0.1) {
	phi = r45;
      } else {
	phi = 0.0;
      }
      cyl.r = r; cyl.phi = phi; cyl.z = z;
      cart = cyl_to_cart(cyl);

      if (vel) {

        if (drift_velocity(cart, 1.0, &v, &setup) == 0) {  // 1.0 => hole charge
          fprintf(file, "%6.5f %6.5f %7.2f\n",
                  r,  z, 100.0*sqrt(v.x*v.x + v.y*v.y + v.z*v.z));
        } else if (pad) {
          fprintf(file, "%6.5f %6.5f %7.2f\n",
                  r,  z, 0.0);
        }

      } else {

        setup.verbosity = TERSE;
        seg = get_signal(cart, s, &setup);
        if (seg < 0) {
          if (0) printf("point not in crystal: (x = %.1f, y = %.1f, z = %.1f)\n", 
                        cart.x, cart.y, cart.z);
          if (pad) fprintf(file, "%6.2f %6.2f %5d\n", r,  z,  0);
          old = -1; // set negative to indicate that we should not yet look for isochrones
        } else {
          n = drift_path_h(&dpath_h, &setup);
          for (j = 0; j < n &&
                 (dpath_h[j].z > setup.pc_length ||
                  (dpath_h[j].x*dpath_h[j].x + dpath_h[j].y*dpath_h[j].y) >
                   setup.pc_radius*setup.pc_radius); j++) ;
          new = (float) j * setup.step_time_calc;
          fprintf(file, "%6.2f %6.2f %5.0f\n", r,  z,  new);

          if (old > 0 && z > 1.1 &&
              (z > setup.ditch_depth+0.5 ||
               r > setup.wrap_around_radius+0.5 ||
               r < -setup.wrap_around_radius-0.5)) {
            for (i=1; i<25; i++) {
              f = 100*i;
              // if (i < 10 && z >60) continue;
              if (((new >= f && old < f) ||
                   (new < f && old >= f)) && ni[i] < 1000) {
                ri[i][ni[i]] = r;
                zi[i][ni[i]++] = z - 0.5*(new-f)/(new-old);
              }
            }
          }
          old = new;
        }

      }
    }
    fprintf(file, "\n");
  }
  fclose(file);

  if (!vel) {
    file = fopen("isochrone.dat", "w");
    for (i=1; i<25; i++) {
      fprintf(file, "## isochrone %d ns :\n##   r (mm)   z (mm)\n", 100*i);
      lneg = lpos = -1;
      if (setup.pc_length < 20 && i < 10) {  // pc_length < 20 to avoid this code for coax dets
        for (j=0; j<ni[i]-1; j++) {
          if (ri[i][j] == ri[i][j+1]) {  // isochrone is double-valued in r
            if (ri[i][j] < 0) {
              if (zi[i][j] > lneg) lneg = zi[i][j];
            } else {
              if (zi[i][j] > lpos) lpos = zi[i][j];
            }
          }
        }
      }
      printf("i = %d  lneg, lpos = %.1f, %.1f\n", i, lneg, lpos);
      if (lneg < 0 && lpos < 0) {
        for (j=0; j<ni[i]; j++) {
          fprintf(file, "%6.2f %6.2f\n", ri[i][j], zi[i][j]);
        }
      } else {
        for (j=ni[i]-1; j>=0; j--) {
          if (ri[i][j] < 0 && zi[i][j] <= lneg)
            fprintf(file, "%6.2f %6.2f\n", ri[i][j], zi[i][j]);
        }
        for (j=0; j<ni[i]; j++) {
          if ((ri[i][j] < 0 && zi[i][j] > lneg) ||
              (ri[i][j] >= 0 && zi[i][j] > lpos))
            fprintf(file, "%6.2f %6.2f\n", ri[i][j], zi[i][j]);
        }
        for (j=ni[i]-1; j>=0; j--) {
          if (ri[i][j] >= 0 && zi[i][j] <= lpos)
            fprintf(file, "%6.2f %6.2f\n", ri[i][j], zi[i][j]);
        }
      }
      fprintf(file, "\n");
    }

  }

  return 0;
}
