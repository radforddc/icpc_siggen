/* signal_calc_util.c -- based on m3d2s.f by I-Yang Lee
 * Karin Lagergren
 *
 * This module includes misc. I/O functions
 */

#include <stdio.h>
#include <strings.h>
#include <string.h>
#include <ctype.h>
#include <stdarg.h>

#include "mjd_siggen.h"

/*this is the default error reporting function.*/
static int default_error_func(const char *format, va_list ap){
  return vfprintf(stderr, format, ap);
}

static int (*output_function)(const char *format, va_list ap) = vprintf;
static int (*error_function)(const char *format, va_list ap) 
     = default_error_func;
static int verbosity_level = NORMAL;

/* tell
   Write stdout using output_function, provided that verb_level is
   above the threshold*/
int tell(int verb_level, const char *format, ...){
  va_list ap;

  if (verbosity_level < verb_level)
    return 0;
  va_start(ap, format);

  output_function(format, ap);
  va_end(ap);
  return 0;
}

/*error
  Error messages are reported using the function pointed to by error_function*/
int error(const char *format, ...){
  va_list ap;
  
  va_start(ap, format);
  error_function(format, ap);
  va_end(ap);
  return 0;
}


/*helper function for reading setup files. 
 returns 0 for success*/
int read_setup_line(FILE *fp, char line[MAX_LINE]){
  char *cp, *cp2;

  *line = '\0';
  for ( ; ; ){
    if (fgets(line, MAX_LINE, fp) == NULL)
      return -1;
    if ((cp = index(line, '#')) != NULL)
      *cp = '\0';
    /*find first non-space character*/
    for (cp = line; *cp != '\0' && isspace(*cp); cp++);
    /*check for empty line*/
    if (strlen(cp) == 0) continue;
    /*removed whitespace @ BOL*/
    if (cp != line){
      for (cp2 = line; *cp != '\0'; cp2++)
	*cp2 = *cp++;
      *cp2 = *cp;//Was BUG
    }
    /*remove whitespace @ EOL*/
    for (cp = line + strlen(line)-1; isspace(*cp); cp--) *cp = '\0';
    return 0;
  }
}


/* set_signal_calc_output
   by default, the verbosity level is set to "normal" and
   messages are written to stdout
   Output can be disabled completely by setting these to NULL. 
   The verbosity level only affects stdout.
*/

int set_signal_calc_output(int verbosity, 
			   int (*out_fn)(const char *format, va_list ap)){

  output_function = out_fn;
  verbosity_level = verbosity;
  return 0;
}
/* set_signal_calc_error_output
   error messages are written to stderr by default. 
   error reporting can be turned off by calling this function with 
   a NULL argument. Can be redirected by supplying an
   alternate output function
*/
int set_signal_calc_error_output(int (*err_fn)(const char *format, va_list ap)){
  error_function = err_fn;
  return 0;
}


