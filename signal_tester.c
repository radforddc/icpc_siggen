/* signal_tester.c
 *  Karin Lagergren
 *
 * a simple program that communicates with the signal 
 * calculation code. Given coordinate => spectrum, print-out of signal.
 * "cart" and "cyl" switch between cartesian and cylindrical coordinates
 * for input
 *
 * Modified Oct 2014 David Radford
 *   to use new config file in common with mjd_fieldgen
 *   to remove all global variables
 *   and to remove code not needed for unsegmented point contact detectors
 *
 * to compile: 
 *  gcc -O3 -Wall -o stester calc_signal.c cyl_point.c detector_geometry.c \
 *                   fields.c point.c read_config.c signal_tester.c -lm -lreadline
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <readline/readline.h>
#include <readline/history.h>
#include <ctype.h>
#include <string.h>

#include "mjd_siggen.h"
#include "calc_signal.h"
#include "cyl_point.h"
#include "detector_geometry.h"
#include "fields.h"

#define PROMPT ">> "
#define MAX_LINE 512
#define MAX_SPE_CHS 16384
#define FREE(p) do{free(p); p = NULL; }while(0)

static int get_cart(char *cmd, struct point *cart, int coord_type);
static int rl_gets(char *line, int MAX_LEN);
static int write_spectrum(float *spec, int nchs, char *name);
static int set_cyl(char *cmd, MJD_Siggen_Setup *setup);
static int set_cart(char *cmd, MJD_Siggen_Setup *setup);
static int save_signal(char *cmd, MJD_Siggen_Setup *setup);
static int print_signal(char *cmd, MJD_Siggen_Setup *setup);
static int print_help(char *cmd, MJD_Siggen_Setup *setup);
static int drift_paths(char *cmd, MJD_Siggen_Setup *setup);
static int set_temp_local(char *cmd, MJD_Siggen_Setup *setup);
static int set_tau(char *cmd, MJD_Siggen_Setup *setup);
static int set_charge_size(char *cmd, MJD_Siggen_Setup *setup);
static int set_diffusion(char *cmd, MJD_Siggen_Setup *setup);
static int set_energy(char *cmd, MJD_Siggen_Setup *setup);
static int set_verbosity(char *cmd, MJD_Siggen_Setup *setup);

static struct{
  char cmd[MAX_LINE];
  int (*f)(char *, MJD_Siggen_Setup *);
  char help_str[MAX_LINE];
}cmds[] = {{"cyl", set_cyl, "cyl ; set input system to cylindrical coords"},
           {"cart", set_cart, "cart ; set input system to cartesian coords"},
           {"sig", save_signal, "sig x y z file.spe or sig r p z file.spe ; save signal"},
	   {"psig", print_signal, "psig x y z or psig r p z ; print signal"},
	   {"dp", drift_paths, "dp fn.dat ; extract charge drift paths to fn.dat"},
	   {"st", set_temp_local, "st %f ; set temperture in K"},
	   {"tau", set_tau, "tau %f ; set preamp integration time in ns"},
	   {"ccs", set_charge_size, "ccs %f ; set charge cloud size in mm"},
	   {"dif", set_diffusion, "dif 0/1 ; set diffusion off/on"},
	   {"ene", set_energy, "ene %f ; set interaction energy in keV"},
	   {"verb", set_verbosity, "verb %i ; set verbosity level [0/1/2]"},
           {"help", print_help, "help ; this output\nquit ; exit program"}};

/* ------------------------------------------ */

int main(int argc, char **argv) {

  MJD_Siggen_Setup setup;

  char  ans[MAX_LINE], *cp;
  int   i, ncmds;

  
  /* read config file and initialize */
  if (argc < 2) {
    printf("Usage: %s <config_file_name>\n", argv[0]);
    return 1;
  }
  if (signal_calc_init(argv[1], &setup) != 0) return 1;
  setup.coord_type = CYL;

  ncmds = sizeof(cmds)/sizeof(cmds[0]);
  while (1) {
    rl_gets(ans, MAX_LINE);
    for (cp = ans; isspace(*cp); cp++);
    if (strlen(cp) == 0) continue;
    if (strncmp(cp, "quit", 4) == 0) return 0;
    for (i = 0; i < ncmds; i++){
      if (strncmp(cp, cmds[i].cmd, strlen(cmds[i].cmd)) == 0){
	cp += strlen(cmds[i].cmd);
	cmds[i].f(cp, &setup);
	break;
      }
    } 
    if (i == ncmds){
      printf("unknown command: %s\n", ans);
    }
  }
  return 0;
}

/*readline wrapper*/
static int rl_gets(char *line, int MAX_LEN){
  static char *line_read = (char *) NULL;

  if (line_read != NULL){
    FREE(line_read);
    line_read = (char *) NULL;
  }

  line_read = readline(PROMPT);
  if (line_read != NULL && *line_read != '\0'){
    add_history(line_read);
  }
  strncpy(line,line_read,MAX_LEN);
  line[MAX_LEN -1] = '\0';
  FREE(line_read);

  return 1;
}

struct spe_header{
  int reclA;            /* 24 */
  unsigned int title[2]; /*don't know why this is uint, but seems to
                           work, so...*/ 
  int dim;
  int a1;               /*  1 */
  int a2;               /*  1 */
  int a3;               /*  1 */
  int reclB;            /* 24 */
};

/*write_spectrum
 *
 * saves the spectrum pointed to by "spec". Returns 0 if unsuccessful,
 * 1 if successful.
 * The name is truncated by removing any trailing .spe (and any
 * subsequent characters...), as well as any leading directory names.
 * If the resulting string is longer than the maximum allowed 8
 * characters, only the first 8 are retained
 */
static int write_spectrum(float *spec, int nchs, char *name){  
  FILE *fp;
  int record_length;
  struct spe_header header;
  char *suffix_start;
  char *fname_start;

  header.reclA = header.reclB = 24; 
  header.title[0] = header.title[1] = 0;
  header.a1 = header.a2 = header.a3 = 1;

  fp = fopen(name,"w");
  if (fp == NULL){
    fprintf(stderr,"Error! Unable to open spectrum file %s \n",name);
    return 0;
  }
  header.dim = nchs;
  if ((suffix_start = strstr(name,".spe")) == NULL){
    suffix_start = name + strlen(name);
  }
  if ((fname_start = rindex(name,'/')) == NULL){
    fname_start = name; 
  }else{
    fname_start++;/*get rid of the '/'*/
  }
  if (suffix_start - fname_start < 8){
    memcpy(header.title,"       ",8);/*blank the title*/
    memcpy(header.title,fname_start,suffix_start - fname_start);
  }else{ 
    memcpy(header.title,suffix_start - 8,8);
  }
  record_length = sizeof(float)*header.dim;

  fwrite(&header, sizeof(header), 1, fp);
  fwrite(&record_length, sizeof(record_length), 1, fp);
  fwrite(spec, sizeof(float), nchs, fp); 
  fwrite(&record_length, sizeof(record_length), 1,fp);
  fclose(fp);

  return 1;
}


static int set_cyl(char *cmd, MJD_Siggen_Setup *setup){
  setup->coord_type = CYL;
  printf("coordinate system: cylindrical\n");
  return 0;
}

static int set_cart(char *cmd, MJD_Siggen_Setup *setup){
  setup->coord_type = CART;
  printf("coordinate system: cartesian\n");
  return 0;
}

/* extracts (cartesian) coordinate from cmd
   if coord_type is CYL, extracts cylindrical coord and converts to cartesian
   on error, returns -1
   on success, returns ending position of coordinates inside cmd
 */
static int get_cart(char *cmd, struct point *cart, int coord_type) {
  float x, y, z;
  struct cyl_pt cyl;
  char  *cp, *cp2;

  if (sscanf(cmd, "%f %f %f", &x, &y, &z) != 3 &&
      sscanf(cmd, "%f, %f, %f", &x, &y, &z) != 3) {
    printf("error parsing coordinate: %s\n", cmd);
    return -1;
  }
  if (coord_type == CYL) {
    cyl.r = x; cyl.phi = y; cyl.z = z;
    *cart = cyl_to_cart(cyl);
  } else {
    cart->x = x; cart->y = y; cart->z = z;
  }
  
  for (cp = cmd; isspace(*cp); cp++);
  for ( ; isdigit(*cp) || *cp == '.'; cp++); //skip r/x coord.
  for ( ; isspace(*cp) || *cp == ','; cp++); //skip whitespace, commas
  for ( ; isdigit(*cp) || *cp == '.'; cp++); //skip phi/y coord.
  for ( ; isspace(*cp) || *cp == ','; cp++); //skip whitespace, commas
  for ( ; isdigit(*cp) || *cp == '.'; cp++); //skip z coord
  for ( ; isspace(*cp) || *cp == ','; cp++); //skip whitespace, commas
  for (cp2 = cp + strlen(cp); isspace(*cp2); cp2--) //remove trailing whitespace
    *cp2 = '\0';

  return (int) (cp - cmd);
}

static int save_signal(char *cmd, MJD_Siggen_Setup *setup){
  float spec[MAX_SPE_CHS], d, a_over_e=0, t10, t90;
  struct point cart;
  int   comp_f;
  int   i, j, k, a_over_e_dt=4;
  char  *cp;
  static float *s;

  if (s == NULL){//first call
    if ((s = (float *) malloc(setup->ntsteps_out*sizeof(*s))) == NULL){
      printf("Malloc failed\n");
      return 1;
    }
  }

  if ((k = get_cart(cmd, &cart, setup->coord_type)) <= 0) return 1;
  cp = cmd + k;
  if (strlen(cp) == 0){
    fprintf(stderr, "must supply file name\n");
    return 1;
  }

  if (get_signal(cart, s, setup) < 0) {
    printf("point not in crystal or has no field: (x = %.1f, y = %.1f, z = %.1f)\n",
	   cart.x, cart.y, cart.z);
    return 1;
  }

  /* calculate t10, t90, and A/E */
  a_over_e = t10 = t90 = 0;
  for (i = 0; i < setup->ntsteps_out - a_over_e_dt; i++) {
    d = (s[i+a_over_e_dt] - s[i]) / (float) a_over_e_dt;
    if (d > a_over_e) a_over_e = d;
    if (s[i] < 0.1f && s[i+1] >= 0.1f)
      t10 = 10.0f*((float) i + (0.1f - s[i]) / (s[i+1] - s[i]));
    if (s[i] < 0.9f && s[i+1] >= 0.9f)
      t90 = 10.0f*((float) i + (0.9f - s[i]) / (s[i+1] - s[i]));
  }
  printf("t_90, t_10-90, A/E = %4.0f %4.0f %f\n", t90, t90-t10, a_over_e);

  comp_f = 1;
  while (setup->ntsteps_out/comp_f > MAX_SPE_CHS) comp_f++;
  if (comp_f > 1) printf("spectrum will be compressed by factor %d\n", comp_f);
  for (i = 0; i < (int) (sizeof(spec)/sizeof(spec[0])); i++) spec[i] = 0.0;

  /*copy signal data to spectrum array*/
  for (j = 0; j < setup->ntsteps_out; j++) spec[j/comp_f] += s[j]*1000/comp_f;
  write_spectrum(spec, setup->ntsteps_out/comp_f, cp);
  printf("%d channels of data saved in spectrum %s\n",
	 setup->ntsteps_out/comp_f, cp);

  return 0;
}

static int print_signal(char *cmd, MJD_Siggen_Setup *setup){
  struct point cart;
  int    i;
  static float *s;

  if (s == NULL){//first call
    if ((s = (float *) malloc(setup->ntsteps_out*sizeof(*s))) == NULL) {
      printf("Malloc failed\n");
      return 1;
    }
  }

  if (get_cart(cmd, &cart, setup->coord_type) <= 0) return 1;

  if (get_signal(cart, s, setup) < 0) {
    printf("point not in crystal or has no field: (x = %.1f, y = %.1f, z = %.1f)\n",
	   cart.x, cart.y, cart.z);
    return 1;
  }

  printf("signal: \n");
  for (i = 0; i < setup->ntsteps_out; i++){
    printf("%.3f ", s[i]);
    if (i%10 == 9) printf("\n");
  }
  printf("\n");
  return 0;
}

static int set_charge_size(char *cmd, MJD_Siggen_Setup *setup){
  float cs;
  char *endp;

  cs = strtod(cmd, &endp);
  if (endp == cmd){
    printf("cannot parse charge_cloud_size: %s\n", cmd);
    return 1;
  }
  if (cs < 0) cs = 0;
  setup->charge_cloud_size = cs;
  printf("Charge cloud size set to %f mm\n", setup->charge_cloud_size);

  return 0;
}

static int set_diffusion(char *cmd, MJD_Siggen_Setup *setup){
  int  d;
  char *endp;

  d = (int) strtol(cmd, &endp, 0);
  if (endp == cmd){
    printf("cannot parse diffusion option: %s\n", cmd);
    return 1;
  }
  if (d < 0) d = 0;
  setup->use_diffusion = d;
  if (d > 0) {
    printf("Diffusion turned on\n");
  } else {
    printf("Diffusion turned off\n");
  }

  return 0;
}

static int set_energy(char *cmd, MJD_Siggen_Setup *setup){
  float e;
  char *endp;

  e = strtod(cmd, &endp);
  if (endp == cmd){
    printf("cannot parse energy: %s\n", cmd);
    return 1;
  }
  if (e < 0) e = 0;
  setup->energy = e;
  printf("Interaction energy set to %.1f keV\n", e);

  return 0;
}

static int set_tau(char *cmd, MJD_Siggen_Setup *setup){
  float t;
  char *endp;

  t = strtod(cmd, &endp);
  if (endp == cmd){
    printf("cannot parse tau: %s\n", cmd);
    return 1;
  }
  if (t < 0) t = 0;
  setup->preamp_tau = t;
  printf("Signals will be integrated with tau = %f ns\n", setup->preamp_tau);

  return 0;
}

static int set_temp_local(char *cmd, MJD_Siggen_Setup *setup){
  float t;
  char *endp;
  
  t = strtod(cmd, &endp);
  if (endp == cmd){
    printf("cannot parse temperature: %s\n", cmd);
    return 1;
  }
  set_temp(t, setup);
  return 0;
}

static int print_help(char *cmd, MJD_Siggen_Setup *setup){
  int i;

  printf("available commands:\n");
  for (i = 0; i < (int) (sizeof(cmds)/sizeof(cmds[0])); i++){
    printf("%s\n",cmds[i].help_str);
  }
  return 0;
}

static int drift_paths(char *cmd, MJD_Siggen_Setup *setup){
  FILE *fp;
  point *dpe, *dph;
  int i, nt;
  char *cp, *cp2;

  int h_only=0;  // set to 1 for use with speed_path_iso figures


  if ((cp = strstr(cmd, "-h"))) {
    while (*(cp-1) == ' ') cp--;
    *cp = 0;
    h_only = 1;
  }

  for (cp = cmd; isspace(*cp); cp++);
  for (cp2 = cp + strlen(cp); isspace(*cp2); cp2--)//remove whitespace
    *cp2 = '\0';
  if (strlen(cp) == 0){
    fprintf(stderr, "must supply file name\n");
    return 1;
  }

  if ((fp = fopen(cp, "w")) == NULL){
    printf("failed to open output file: %s\n", cp);
    return 1;
  }

  if (h_only) {
    nt =  drift_path_h(&dph, setup);
    fprintf(fp,"#hx hy hz\n");
    for (i = 0; i < nt && dph[i].z > 0.0; i++){
      fprintf(fp,"%7.3f %7.3f %7.3f\n", dph[i].x, dph[i].z, -dph[i].x);
    }
  } else {
    nt = drift_path_e(&dpe, setup);
    drift_path_h(&dph, setup);
    fprintf(fp,"# n ex ey ez hx hy hz\n");
    for (i = 0; i < nt; i++) {
      if (dpe[i].x == 0 && dpe[i].y == 0 && dpe[i].y == 0 &&
          dph[i].x == 0 && dph[i].y == 0 && dph[i].y == 0) break;
      fprintf(fp,"%4d %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f\n",
              i, dpe[i].x, dpe[i].y, dpe[i].z, dph[i].x, dph[i].y, dph[i].z);
    }
  }

  fprintf(fp,"\n");
  fclose(fp);
  return 0;
}

static int set_verbosity(char *cmd, MJD_Siggen_Setup *setup) {
  int i;
  char *endp;

  i = strtol(cmd, &endp, 0);
  if (endp == cmd) {
    printf("cannot parse verbosity level: %s\n", cmd);
    return 1;
  }
  if (i < 0) i = 0;
  if (i > 2) i = 2;
  setup->verbosity = i;
  printf("Verbosity level set to %d\n", i);

  return 0;
}
