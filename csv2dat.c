#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

struct spe_header{
  int reclA;               /* 24 */
  char title[8];
  int dim;
  int a1, a2, a3, reclB;   /*  1 1 1 24 */
};

static int write_spectrum(float *spec, int nchs, char *name) {  
  FILE *fp;
  int record_length;
  struct spe_header header;

  header.reclA = header.reclB = 24; 
  header.title[0] = header.title[1] = 0;
  header.a1 = header.a2 = header.a3 = 1;

  fp = fopen(name,"w");
  if (fp == NULL){
    fprintf(stderr,"Error! Unable to open spectrum file %s \n",name);
    return 1;
  }
  header.dim = nchs;
  memcpy(header.title, name, 8);
  record_length = sizeof(float)*header.dim;

  fwrite(&header, sizeof(header), 1, fp);
  fwrite(&record_length, sizeof(record_length), 1, fp);
  fwrite(spec, sizeof(float), nchs, fp); 
  fwrite(&record_length, sizeof(record_length), 1,fp);
  fclose(fp);

  return 0;
}

int main(int argc, char **argv) {

  int    i, nchs = 0, column = 0, reverse = 0;
  char   line[120];
  float  spec[16384] = {0}, spec2[16384], f[5], factor = -1;  // NOTE default factor is -1, for p-type detectors in mjd_fieldgen
  FILE  *fp;

  if (argc < 4) {
    printf("Usage: %s <csv_file_name> <column_number> <output_file_name> [<optional_factor>] [-r]\n"
           "       Use .spe or .dat as the output file name extension\n"
           "       Column number starts at 1\n"
           "       Default factor value is -1, for p-type detectors in mjd_fieldgen\n"
           "       -r: Reverse ordering of spectrum for use with mjd_fieldgen\n", argv[0]);
    return 1;
  }
  if (argc > 4) factor = atof(argv[4]);
  if (factor == 0) factor = -1;
  if ((argc > 4 && strstr(argv[4], "-r")) ||
      (argc > 5 && strstr(argv[5], "-r"))) reverse = 1;
  
  if ((fp = fopen(argv[1], "r")) == NULL) {
    printf("Error: Cannot open input file %s\n", argv[1]);
    return 1;
  }
  // skip first line of file (one-line header)
  // fgets(line, 120, fp);

  column = atoi(argv[2]);
  if (column < 1 || column> 5) {
    printf("Error: Invalid column number %d\n", column);
    return 1;
  }

  /* read data */
  printf("\nReading file %s, column %d, factor %f, reverse = %d\n", argv[1], column, factor, reverse);
  while (nchs < 16384 && fgets(line, 120, fp)) {
    if (sscanf(line, "%f,%f,%f,%f,%f", f, f+1, f+2, f+3, f+4) >= column) {
      spec[nchs++] = f[column - 1] * factor;
    }
  }
  fclose(fp);

  printf("%i channels read\n", nchs);
  if (nchs < 4) {
    printf("Cannot read input file %s\n", argv[1]);
    return 1;
  }
  /* ------------  if profile slope has wrong direction,
                   then reverse ordering of spectrum for use with mjd_fieldgen ---------- */
  if ( (reverse || fabs(spec[0]) < fabs(spec[nchs-1])) &&
      !(reverse && fabs(spec[0]) < fabs(spec[nchs-1]))) {
    memcpy(spec2, spec, sizeof(spec));
    for (i=0; i<nchs; i++) spec[i] = spec2[nchs - i - 1];
  }

  if (nchs < 200) nchs = 200;

  /* write .spe or .dat file */
  if (strstr(argv[3], ".spe")) {
    if (write_spectrum(spec, nchs, argv[3])) return 1;
  } else if (strstr(argv[3], ".dat")) {
    if (!(fp = fopen(argv[3], "w"))) {
      printf("Error: Cannot open ioutput file %s!\n", argv[3]);
      return 1;
    }
    fwrite(spec, sizeof(float), nchs, fp); 
    fclose(fp);
  } else {
    printf("Error: Invalid file name extension in %s!\n", argv[3]);
    return 1;
  }

  printf("\n Success. Wrote %d chs to %s\n\n", nchs, argv[3]);
  return 0;
}
