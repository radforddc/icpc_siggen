/* drift_time_map.c
 *  D.C. Radford
 *
 * code to generate a map of drift times (in ns) for point-contact detectors
 * uses t90 as drift time, assuming that siggen output time steps are 10 ns
 *
 * to compile: 
 *  gcc -o drift_time_map drift_time_map.c calc_signal.c cyl_point.c \
 *  detector_geometry.c fields.c point.c read_config.c drift_time_map.c -lm
 *
 * to run: drift_time_map <config_file> <DT_map_output_file>
 *
 * it might be good to set verbosity_level to 0 in the config file,
 * and to set tau to 0 and step_time_out to 10.0
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

/* ------------------------------------------ */

int main(int argc, char **argv) {

  MJD_Siggen_Setup setup;
  int   i, r, z;
  FILE *fp;
  struct point cart;
  float *s, t90;

  
  /* read config file and initialize */
  if (argc < 3) {
    printf("Usage: %s <config_file_name> <DT_map_output_filename>\n", argv[0]);
    return 1;
  }
  if (signal_calc_init(argv[1], &setup) != 0) return 1;
  setup.coord_type = CYL;

  /* open output file */
  if (!(fp=fopen(argv[2], "w"))) {
    printf("Error, cannot open output file %s\n", argv[2]);
    return 1;
  }

  /* malloc space for signal */
  if ((s = (float *) malloc(setup.ntsteps_out*sizeof(*s))) == NULL) {
    printf("Malloc failed\n");
    return 1;
  }
  fprintf(fp, "#  r   z    drift_time (t90)\n");

  /* step over r and z */
  for (r=0; r < lrint(setup.xtal_radius); r++) { //use one-mm steps in radius
    for (z=0; z < lrint(setup.xtal_length); z++) { //use one-mm steps in length
      cart.x = r; cart.y = 0; cart.z = z;
      t90 = 0;
      if (get_signal(cart, s, &setup) >= 0) {
        for (i = 0; i < setup.ntsteps_out-1; i++) {
          if (s[i] < 0.9f && s[i+1] >= 0.9f) break;
        }
        t90 = 10.0f*((float) i);
        if (i < setup.ntsteps_out-1) {
          t90 += 10.0f*((0.9f - s[i]) / (s[i+1] - s[i]));
        }
      }
      fprintf(fp, " %3d %3d %6.1f\n", r, z, t90);
    }
  }

  fclose(fp);
  return 0;
}
