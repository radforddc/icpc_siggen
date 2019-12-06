/*
 *  mjd_siggen.h -- define data structures used by rewritten fieldgen and siggen for MJD detectors
 *                    (BEGes and PPC detectors)
 *  David Radford   Oct 2014
 */

#ifndef _MJD_SIGGEN_H
#define _MJD_SIGGEN_H

#include "point.h"
#include "cyl_point.h"

/* verbosity levels for std output */
#define TERSE  0
#define NORMAL 1
#define CHATTY 2

// #define VERBOSE 2  // Set to 0 for quiet, 1 or 2 for less or more info
#define TELL_NORMAL if (setup->verbosity >= NORMAL) tell
#define TELL_CHATTY if (setup->verbosity >= CHATTY) tell

/* Reference temperature for drift vel. corrections is 77K */
#define REF_TEMP 77.0
/* max, min temperatures for allowed range */
#define MIN_TEMP 40.0
#define MAX_TEMP 120.0
/* enum to identify cylindrical or cartesian coords */
#define CYL 0
#define CART 1
#define MAX_FNAME_LEN 512

float sqrtf(float x);
float fminf(float x, float y);

// from fields.c
struct velocity_lookup{
  float e;
  float e100;
  float e110;
  float e111;
  float h100;
  float h110;
  float h111;
  float ea; //coefficients for anisotropic drift 
  float eb;
  float ec;
  float ebp;
  float ecp;
  float ha;
  float hb;
  float hc;
  float hbp;
  float hcp;
  float hcorr;
  float ecorr;
};

/* setup parameters data structure */
typedef struct {
  // general
  int verbosity;              // 0 = terse, 1 = normal, 2 = chatty/verbose

  // geometry
  float xtal_length;          // z length
  float xtal_radius;          // radius
  float pc_length;            // point contact length
  float pc_radius;            // point contact radius
  float wrap_around_radius;   // wrap-around radius for BEGes. Set to zero for ORTEC
  float ditch_depth;          // depth of ditch next to wrap-around for BEGes. Set to zero for ORTEC
  float ditch_thickness;      // width of ditch next to wrap-around for BEGes. Set to zero for ORTEC
  float bottom_taper_length;  // size of 45-degree taper at bottom of ORTEC-type crystal
  float hole_length;          // length of hole, for inverted-coax style
  float hole_gap;             // length of xtal - length of hole, for inverted-coax style
  float hole_radius;          // radius of hole, for inverted-coax style
  float outer_taper_length;   // z-length of outside taper for inverted-coax style
  float taper_angle;          // taper angle in degrees, for inner or outer taper
  float inner_taper_length;   // z-length of inside (hole) taper for inverted-coax style
  float outer_taper_width;    // r-width of outside taper at far top of crystal
  float inner_taper_width;    // r-width of inside (hole) taper at far top of crystal
  float top_bullet_radius;    // bulletization radius at top of crystal
  float bottom_bullet_radius; // bulletization radius at bottom of BEGe crystal
  float hole_bullet_radius;   // bulletization radius at bottom of hole
  float Li_thickness;         // depth of full-charge-collection boundary for Li contact
  float pc_offset_x;          // x-offset of point-contact position
  float pc_offset_y;          // y-offset of point-contact position
  float pc_offset_max;        // maximum radial offset of pc position
  float hole_offset_x_top;       // x-offset of well position at top of hole
  float hole_offset_y_top;       // y-offset of well position at top of hole
  float hole_offset_x_bottom;    // x-offset of well position at bottom of hole
  float hole_offset_y_bottom;    // y-offset of well position at bottom of hole
  float hole_offset_max;         // maximum radial offset of well position

  // electric fields & weighing potentials
  float xtal_grid;            // grid size in mm for field files (either 0.5 or 0.1 mm)
  float impurity_z0;          // net impurity concentration at Z=0, in 1e10 e/cm3
  float impurity_gradient;    // net impurity gradient, in 1e10 e/cm4
  float impurity_quadratic;   // net impurity difference from linear, at z=L/2, in 1e10 e/cm3
  float impurity_surface;     // surface impurity of passivation layer, in 1e10 e/cm2
  float impurity_radial_add;  // additive radial impurity at outside radius, in 1e10 e/cm3
  float impurity_radial_mult; // multiplicative radial impurity at outside radius (neutral=1.0)
  float impurity_rpower;      // power for radial impurity increase with radius
  float xtal_HV;              // detector bias for fieldgen, in Volts
  int   max_iterations;       // maximum number of iterations to use in mjd_fieldgen
  int   write_field;          // set to 1 to write V and E to output file, 0 otherwise
  int   write_WP;             // set to 1 to calculate WP and write it to output file, 0 otherwise
  int   bulletize_PC;         // set to 1 for inside of point contact hemispherical, 0 for cylindrical
  int   fix_adaptive;         // 

  // file names
  char drift_name[256];       // drift velocity lookup table
  char field_name[256];       // potential/efield file name
  char wp_name[256];          // weighting potential file name
  char config_file_name[256];

  // signal calculation 
  float xtal_temp;            // crystal temperature in Kelvin
  float preamp_tau;           // integration time constant for preamplifier, in ns
  int   time_steps_calc;      // number of time steps used in calculations
  float step_time_calc;       // length of time step used for calculation, in ns
  float step_time_out;        // length of time step for output signal, in ns
  //    nonzero values in the next few lines significantly slow down the code
  float charge_cloud_size;    // initial FWHM of charge cloud, in mm; set to zero for point charges
  int   use_diffusion;        // set to 0/1 for ignore/add diffusion as the charges drift
  float energy;               // set to energy > 0 to use charge cloud self-repulsion, in keV

  int   coord_type;           // set to CART or CYL for input point coordinate system
  int   ntsteps_out;          // number of time steps in output signal

  // data for fields.c
  float xmin, xmax;
  float ymin, ymax;
  float zmin, zmax;
  int   numx, numy, numz;     // dimensions of efld and wpot arrays

  //for fieldgen:
  double *impurity_z;
  double ***v[2];
  char   ***point_type;
  double ***epsilon;
  float  impurity_lamda, rho_b, rho_c;

  //for siggen:
  int   v_lookup_len;
  struct velocity_lookup *v_lookup;
  float  ***e;
  vector ***efld;
  double ***wpot;

  
  // data for calc_signal.c
  point *dpath_e, *dpath_h;        // electron and hole drift paths
  float surface_drift_vel_factor;  // ratio of velocity on passivated surface rather than in bulk
  float initial_vel, final_vel;    // initial and final drift velocities for charges collected to PC
  float dv_dE;                     // derivative of drift velocity with field ((mm/ns) / (V/cm))
  float v_over_E;                  // ratio of drift velocity to field ((mm/ns) / (V/cm))
  double final_charge_size;        // in mm

} MJD_Siggen_Setup;


int read_config(char *config_file_name, MJD_Siggen_Setup *setup);

// defined in field_init.c:

int geometry_init(MJD_Siggen_Setup *setup);     /* returns # contacts */
enum point_types{OUTSIDE, CONTACT_0, CONTACT_VB, INSIDE, FIXED, PASSIVE, DITCH};
int init_ev_calc(MJD_Siggen_Setup *setup);
int init_wp_calc(MJD_Siggen_Setup *setup);
/* set crystal type. Can be either CRYSTAL_A or CRYSTAL_B
   invalid type => nothing happens. Returns new value
*/
int set_crystal_geometry(MJD_Siggen_Setup *setup);

/*  defined in signal_calc_util.c */

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


#endif /*#ifndef _MJD_SIGGEN_H */
