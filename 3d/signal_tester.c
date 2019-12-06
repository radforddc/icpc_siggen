/* signal_tester.c
 *  Karin Lagergren
 *
 * a simple program that communicates with the signal 
 * calculation code. Given coordinate => spectrum, print-out of signal,
 * or segment number. Given two sets of coords=>rms distance
 * "cart" and "cyl" switch between cartesian and cylindrical coordinates
 * for input
 *
 * to compile: 
 *  gcc -o st signal_tester.c point.c cyl_point.c calc_signal.c\
 *    fields.c detector_geometry.c signal_calc_util.c -lm -lreadline
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
 double trunc(double x);
#include <readline/readline.h>
#include <readline/history.h>
#include <ctype.h>

#include "mjd_siggen.h"
#include "calc_signal.h"
#include "cyl_point.h"
#include "signal_calc_util.h"
#include "detector_geometry.h"
#include "fields.h"

#define PROMPT ">> "
#define RZ_STEP 0.5
#define XYZ_STEP 0.5
#define P_STEP (M_PI/100)
#define MAX_LINE 512
#define MAX_SPE_CHS 16384

#define FREE(p) do{free(p); p = NULL; }while(0)

#define CYL 0
#define CART 1

#define C2_SIGNAL_MULT_F 10000
#define C2_THRESH 100
#define NON_LINEAR_CHISQ 1

static int coord_type = CART;

static int rl_gets(char *line, int MAX_LEN);
static int write_spectrum(float *spec, int nchs, char *name);
// static int rc_integrate(float **s_in, float **s_out, float tau);
int        shift_signal(float **s, float dt, int nsigs, int ntimes);

static int save_signal(MJD_Siggen_Setup *setup, char *cmd);
static int print_help(MJD_Siggen_Setup *setup, char *cmd);
static int aoe(MJD_Siggen_Setup *setup, char *cmd);
static int drift_paths(MJD_Siggen_Setup *setupp, char *cmd);

static int time_steps;

static struct{
  char cmd[MAX_LINE];
  int (*f)(MJD_Siggen_Setup *, char *);
  char help_str[MAX_LINE];
}cmds[] = {{"sig", save_signal, "sig x y z file.spe or sig r p z file.spe"},
           {"aoe", aoe, "aoe r dp z"},
	   {"dp", drift_paths, "dp fn.dat ; extract charge drift paths to fn.dat"},
	   {"help", print_help, "help\nquit ; exit program"}};

int main(int argc, char **argv)
{

  MJD_Siggen_Setup setup;
  char ans[MAX_LINE], *cp;
  int i, ncmds;
  
#ifdef CHATTY
  set_signal_calc_output(CHATTY, vprintf);
#endif

  /* read config file and initialize */
  if (argc < 2) {
    printf("Usage: %s <config_file_name>\n", argv[0]);
    return 1;
  }
  strncpy(setup.config_file_name, argv[1], 256);
  if (signal_calc_init(&setup) != 0) return 1;
  time_steps = setup.ntsteps_out;

  ncmds = sizeof(cmds)/sizeof(cmds[0]);
  while (1) {
    rl_gets(ans, MAX_LINE);
    for (cp = ans; isspace(*cp); cp++);
    if (strlen(cp) == 0) continue;
    if (strstr(ans, "quit")) return 0;
    for (i = 0; i < ncmds; i++){
      if (strncmp(cp, cmds[i].cmd, strlen(cmds[i].cmd)) == 0){
	cp += strlen(cmds[i].cmd);
	cmds[i].f(&setup, cp);
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

static int save_signal(MJD_Siggen_Setup *setup, char *cmd){
  float spec[MAX_SPE_CHS];
  struct point cart;
  struct cyl_pt cyl;
  int comp_f;
  int i, j;
  char *cp, *cp2;
  static float *s;

  if (s == NULL) {//first call
    if ((s = malloc(time_steps*sizeof(*s))) == NULL){
      printf("malloc failed\n");
      return 1;
    }
  }

  if (coord_type == CYL ){
    if (sscanf(cmd, "%f %f %f", &cyl.r, &cyl.phi, &cyl.z) == 3
	  || sscanf(cmd, "%f, %f, %f", &cyl.r, &cyl.phi, &cyl.z) == 3){
      cart = cyl_to_cart(cyl);
    } else {
      printf("error parsing coordinate: %s\n", cmd);
      return 1;
    }
  } else {
    if (sscanf(cmd, "%f %f %f", &cart.x, &cart.y, &cart.z) == 3
	  || sscanf(cmd, "%f, %f, %f", &cart.x, &cart.y, &cart.z) == 3){
      ;//nothing
    } else {
      printf("error parsing coordinate: %s\n", cmd);
      return 1;
    }
  }
  
  for (cp = cmd; isspace(*cp); cp++);
  for ( ; isdigit(*cp) || *cp == '.' || *cp == '-'; cp++); //skip r/x coord.
  for ( ; isspace(*cp) || *cp == ',' || *cp == '-'; cp++); //skip whitespace, commas
  for ( ; isdigit(*cp) || *cp == '.' || *cp == '-'; cp++); //skip phi/y coord.
  for ( ; isspace(*cp) || *cp == ',' || *cp == '-'; cp++); //skip whitespace, commas
  for ( ; isdigit(*cp) || *cp == '.' || *cp == '-'; cp++); //skip z coord

  for ( ; isspace(*cp); cp++); //skip whitespace
  for (cp2 = cp + strlen(cp); isspace(*cp2); cp2--)//remove whitespace
    *cp2 = '\0';
  if (strlen(cp) == 0) {
    fprintf(stderr, "must supply file name\n");
    return 1;
  }

  if (get_signal(setup, cart, s) < 0) {
    printf("point not in crystal: (x = %.1f, y = %.1f, z = %.1f)\n", 
	   cart.x, cart.y, cart.z);
    return 1;
  }
  for (comp_f = 1; time_steps/comp_f > MAX_SPE_CHS; comp_f *= 2)
    ;

  printf("spectrum will be compressed by factor %d\n", comp_f);
  for (i = 0; i < sizeof(spec)/sizeof(spec[0]); i++)
    spec[i] = 0.0;
  /*copy signal data to spectrum array*/
  for (j = 0; j < time_steps; j++){
    spec[j/comp_f] += s[j]*1000/comp_f;
  }
  write_spectrum(spec, time_steps/comp_f, cp);
  printf("%d channels of data saved in spectrum %s\n", time_steps/comp_f, cp);
  float s1 = 0, s2 = 0;
  for (i=0; i < 40; i++) s1 += spec[i+40] - spec[i];
  for (i=0; i < time_steps/comp_f - 80; i++) {
    s1 += spec[i+80] - spec[i+40] - spec[i+39] + spec[i];
    if (s2 < s1) s2 = s1;
  }
  printf(" +><+><+ %.3f %.1f\n", cyl.phi, s2/40.0);
  
  return 0;
}

static int aoe(MJD_Siggen_Setup *setup, char *cmd){
  struct point cart;
  struct cyl_pt cyl;
  int i;
  float dphi, s1, s2;
  static float *s;
  FILE *fp;

  if (s == NULL) {//first call
    if ((s = malloc(time_steps*sizeof(*s))) == NULL){
      printf("malloc failed\n");
      return 1;
    }
  }

  coord_type = CYL;
  if (sscanf(cmd, "%f %f %f",   &cyl.r, &dphi, &cyl.z) != 3 &&
      sscanf(cmd, "%f, %f, %f", &cyl.r, &dphi, &cyl.z) != 3) {
    printf("error parsing coordinate: %s\n", cmd);
    return 1;
  }

  fp = fopen("aoe.dat", "w");
  fprintf(fp, "# r = %.1f, z = %.1f\n#  phi    A/E     DT\n", cyl.r, cyl.z);
  for (cyl.phi = 0; cyl.phi < 6.28; cyl.phi += dphi) {
    cart = cyl_to_cart(cyl);
    if (get_signal(setup, cart, s) < 0) {
      printf("point not in crystal: (x = %.1f, y = %.1f, z = %.1f)\n", 
             cart.x, cart.y, cart.z);
      return 1;
    }

    s1 = s2 = 0;
    //    for (i=0; i < 40; i++) s1 += s[i+40] - s[i];
    //    for (i=0; i < time_steps - 80; i++) {
    //      s1 += s[i+80] - 2.0*s[i+40] + s[i];
    for (i=0; i < 20; i++) s1 += s[i+20] - s[i];
    for (i=0; i < time_steps - 40; i++) {
      s1 += s[i+40] - 2.0*s[i+20] + s[i];
      if (s2 < s1) s2 = s1;
    }
    for (i=0; i < time_steps; i++) {
      if (s[i] > 0.7) break;
    }
    fprintf(fp, " %6.3f %6.1f  %4d\n", cyl.phi, 3.1*s2*1000.0/40.0, i);
  }
  fclose(fp);

  return 0;
}

static int print_help(MJD_Siggen_Setup *setup, char *cmd){
  int i;

  printf("available commands:\n");
  for (i = 0; i < sizeof(cmds)/sizeof(cmds[0]); i++){
    printf("%s\n",cmds[i].help_str);
  }
  return 0;
}


static int drift_paths(MJD_Siggen_Setup *setup, char *cmd){
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
    nt =  drift_path_h(setup, &dph);
    fprintf(fp,"#hx hy hz\n");
    for (i = 0; i < nt && dph[i].z > 0.0; i++){
      fprintf(fp,"%7.3f %7.3f %7.3f\n", dph[i].x, dph[i].z, -dph[i].x);
    }
  } else {
    nt = drift_path_e(setup, &dpe);
    drift_path_h(setup, &dph);
    fprintf(fp,"# n ex ey ez hx hy hz\n");
    for (i = 0; i < nt && (dph[i].z > 0.0 || dpe[i].z > 0.0); i++) {
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
