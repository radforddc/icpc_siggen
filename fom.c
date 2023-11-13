
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


int main(int argc, char **argv) {
  int   i, j;
  int   fom[6][8192] = {{0}};
  FILE  *file_out;

  file_out = fopen("fom.sec", "w");

  for (j=0; j<6; j++) {
    for (i=2000; i<8000; i++)
      fom[j][i] = 1000;
    for (i=0; i<3999-200*j; i++)
      fom[j][i] = 200.0 + 10.0*j + 2000.0/(4000.0 - (float) i - 200.0*j);
  }

  fwrite(fom, sizeof(fom), 1, file_out);
  
  fclose (file_out);

  return 0;
}


