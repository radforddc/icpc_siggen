/* calc_signal.h -- based on m3d2s.f by I-Yang Lee
 * Karin Lagergren
 *
 * This module contains the main interface to the signal calculation
 * code. 
 *
 * To use: 
 * -- call signal_calc_init. This will initialize geometry, fields,
 *       drift velocities etc.
 * -- call hit_segment/get_signal
 */
#ifndef _CALC_SIGNAL_H
#define _CALC_SIGNAL_H

#include <stdarg.h>
#include "mjd_siggen.h"
#include "point.h"

#define MAX_LINE            512
//#define NET_SIGNAL_THRESH   0.55
//#define WP_THRESH           0.55
//#define WP_THRESH_ELECTRONS 1e-4 /* electrons are considered collected if */
                                   /* they stop drifting where the wp is < this */

typedef struct{
  float **s;
  int *t_lo;
  int *t_hi;
} Signal;

/* signal_calc_init
   read setup from configuration file,
   then read the electric field and weighting potentials,
   and initialize the signal calculation variables and nsegments
   returns 0 for success
*/
int signal_calc_init(MJD_Siggen_Setup *setup);

/* hit_segment
 * return the segment number for an interaction at point pt
 * -1 if outside crystal
 */
int hit_segment(MJD_Siggen_Setup *setup, point pt);

/* get_signal calculate signal for point pt. Result is placed in
 * signal array which is assumed to have at least (number of segments)
 * x ( number of time steps) elements
 * returns segment number or -1 if outside crystal
 */
int get_signal(MJD_Siggen_Setup *setup, point pt, float *signal);

/* signal_calc_finalize
 * Clean up
 */
int signal_calc_finalize(MJD_Siggen_Setup *setup);

/*drift paths for last calculated signal.
  after the call, "path" will point at a 1D array containing the points
  (one per time step) of the drift path. 
  freeing that pointer will break the code.
*/
int drift_path_e(MJD_Siggen_Setup *setup, point **path);
int drift_path_h(MJD_Siggen_Setup *setup, point **path);

#endif /*#ifndef _CALC_SIGNAL_H*/
