/* signal_calc_util.h -- based on m3d2s.f by I-Yang Lee
 * Karin Lagergren
 *
 * This module includes misc. I/O functions
 */

#ifndef _UTIL_H
#define _UTIL_H

#include <stdarg.h>
#define MAX_LINE 512

/* read_setup_line
 * read one line from config file
 * # (lumberyard/pound sign/hash) turns rest of line into comment
 * empty lines are skipped, whitespace stripped at beginning and end of line
 * returns 0 for success
 */
int read_setup_line(FILE *fp, char line[MAX_LINE]);

/* set_signal_calc_output
   by default, the verbosity level is set to "normal" and
   messages are written to stdout
   Output can be disabled completely by setting these to NULL. 
   The verbosity level only affects stdout.
   usage example: set_signal_calc_output(TERSE, vprintf);
   any function that matches the given prototype (same as vprintf)
   will work, though.
*/
int set_signal_calc_output(int verbosity, 
			   int (*out_fn)(const char *format, va_list ap));
/* set_signal_calc_error_output
   error messages are written to stderr by default. 
   error reporting can be turned off by calling this function with 
   a NULL argument. Can be redirected by supplying an
   alternate output function
*/
int set_signal_calc_error_output(int (*out_fn)(const char *format, va_list ap));

/* These are the actual functions that are used to print to stdout and stderr,
   respectively. They use the settings above to determine what gets printed
   where
*/
int tell(int verb_level, const char *format, ...);
int error(const char *format, ...);

#endif /*#ifndef _UTIL_H*/
