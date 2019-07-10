/* calc_signal.h -- based on m3d2s.f by I-Yang Lee
 * Karin Lagergren
 *
 * This module contains the main interface to the signal calculation
 * code. 
 *
 * To use: 
 * -- call signal_calc_init. This will initialize geometry, fields,
 *       drift velocities etc.
 * -- call get_signal
 */
#ifndef _CALC_SIGNAL_H
#define _CALC_SIGNAL_H

#include <stdarg.h>
#include "point.h"
#include "mjd_siggen.h"

#define MAX_LINE 512
#define NET_SIGNAL_THRESH 0.55
#define WP_THRESH 0.55
#define WP_THRESH_ELECTRONS 1e-4 /*electrons are considered collected if
				   they stop drifting where the wp is < this*/

typedef struct {
  float *s;
  int   *t_lo;
  int   *t_hi;
} Signal;

/* signal_calc_init
   read setup from configuration file,
   then read the electric field and weighting potential,
   and initialize the signal calculation variables
   returns 0 for success
*/
int signal_calc_init(char *config_file_name, MJD_Siggen_Setup *setup);

/* get_signal calculate signal for point pt. Result is placed in signal
 * array which is assumed to have at least (number of time steps) elements
 * returns -1 if outside crystal
 */
int get_signal(point pt, float *signal, MJD_Siggen_Setup *setup);

/* make_signal
   Generates the signal originating at point pt, for charge q
   returns 0 for success
*/
int make_signal(point pt, float *signal, float q, MJD_Siggen_Setup *setup);

/* signal_calc_finalize
 * Clean up
 */
int signal_calc_finalize(MJD_Siggen_Setup *setup);

/* rc_integrate
 * do RC integratation of signal s_in with time constant tau 
 */
int rc_integrate(float *s_in, float *s_out, float tau, int time_steps);

/*drift paths for last calculated signal.
  after the call, "path" will point at a 1D array containing the points
  (one per time step) of the drift path. 
  freeing that pointer will break the code.
*/
int drift_path_e(point **path, MJD_Siggen_Setup *setup);
int drift_path_h(point **path, MJD_Siggen_Setup *setup);

/* these functions are used to print to stdout and stderr, respectively.
*/
void tell(const char *format, ...);
void error(const char *format, ...);

#endif /*#ifndef _CALC_SIGNAL_H*/
