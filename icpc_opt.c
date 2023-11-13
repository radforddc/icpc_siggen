
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


int main(int argc, char **argv) {
  int   i, j, ng=0, g[5], do_dtc = 0;
  float f, f0, flo[8], fhi[8], fincr[8], fnow[8];
  float x, y, A, B, C, fom, max_fom = -9999, max_f[8], fomm, fomt;
  float wpfrac, fwhm1, fwhm2, aoe1, aoe2;
  float Vlow=5000, Vhigh=0, Elow=400, Ehigh=120, FOMlow=20, FOMhigh=-1;
  float Masslow=5000, Masshigh=1000, FWHMlow=10, FWHMhigh=2, WPlow=100, WPhigh=0, AoElow=10, AoEhigh=0;
  char  *c, line[256], command[256], syscommand[512], emin_loc[32];
  char  title[8][8] = {"wrap_r", "hole_g", "hole_r", "tapr_l", "tapr_a", "xtal_l", "xtal_r", " z_cut"};
  char  ptitle[256];
  FILE  *file_in, *file_out;

  float massgoal=3300, dmass = 200;
  float fwhmgoal= 2.5, dfwhm = 0.15;
  float aoegoal = 0.8, daoe  = 0.3;

  /* read config file and initialize */
  if (argc < 2) {
    printf("Usage: %s <opt_config_file_name> <optional fieldgen flags>\n"
           "          -d  : include DTC calculation (longer run time)\n"
           "          -dd : include diffusion in DTC calculation\n", argv[0]);
    return 1;
  }
  if (!(file_in = fopen(argv[1], "r"))) {
    printf("\nERROR: opt config file %s does not exist?\n", argv[1]);
    return 1;
  }
  for (i=2; i<argc; i++) {
    if (strstr(argv[i], "-d")) {
      do_dtc = 1; // run dtc for drift time correction estmate (slow)
      if (strstr(argv[i], "-dd")) do_dtc = 2; // include charge cloud size, diffusion, self-repulsion
      printf(" Adding DTC calculation...\n");
      argc--;
      for (j=i; j<argc; j++)
        argv[j] = argv[j+1]; // remove -d from argument list
    }
  }
  
  // read root of command to run
  // usually something like "./icpc_fieldgen_opt config_files/mirion_enr/30925.config"
  strncpy(command, "none", 256);
  while (fgets(line, 256, file_in)) {
    if (line[0] == '#') continue;  // comment
    strncpy(command, line, 256);
    break;
  }
  if ((c = strstr(command, "\n"))) *c = 0;
  printf("argc: %d\n",argc);
  for (i=2; i<argc; i++) {
    strcat(command, " ");
    strcat(command, argv[i]);
    printf("%d %s\n",i, argv[i]);
  }
  printf("\n >>> Command root: %s\n", command);
  if (strstr(command, "none")) return 1;
  if ((c = strstr(command, "\n"))) *c = 0;
  // add config file name to title for output plot
  c = strstr(command, " ") + 1;  // start of config file name
  strncpy(ptitle, c, 256);
  *(strstr(ptitle, " ")) = 0;    // end of config file name

  // read geometry options to scan over
  while (fgets(line, 256, file_in) && ng < 5) {
    if (line[0] == '#') continue;  // comment
    if (line[0] == '-' && line[1] == 'd') do_dtc = 1; // run dtc for drift time correction estmate (slow)
    if (line[0] == '-' && line[1] == 'd' && line[2] == 'd') do_dtc = 2; // include charge cloud size, diffusion, self-repulsion
    if (line[0] == '-' && line[1] == 'g') {
      if (*(line+2) >= 'a' && *(line+2) <= 'z') {
        i = -1;
        if (*(line+2) == 'w') i = 0; //   wrap_around_radius
        if (*(line+2) == 'g') i = 1; //   hole_length_gap
        if (*(line+2) == 'h') i = 2; //   hole_radius
        if (*(line+2) == 't') i = 3; //   inner_taper_length
        if (*(line+2) == 'a') i = 4; //   taper_angle
        if (*(line+2) == 'l') i = 5; //   xtal length
        if (*(line+2) == 'r') i = 6; //   xtal radius
        if (*(line+2) == 'z') i = 7; //   z-cut position in mm (adds to any other -z input)
      } else {
        i = atoi(line+2);
      }
      if (i < 0 || i > 7) {
        printf("\nERROR: illegal geometery parameter %s\n", line);
        return -1;
      }
      g[ng] = i;
      if (3 != sscanf(line+3, "%f %f %f", &flo[ng], &fhi[ng], &fincr[ng])) {
        printf("\nERROR: cannot decode geometery parameter list %s\n", line);
        return -1;
      }
      ng++;
      // add flag to title for output plot
      strcat(ptitle, ";  ");
      strcat(ptitle, title[i]);
      if ((c=strstr(line, "#"))) *c=0; // remove trailing comment
      if ((c=strstr(line, "\n"))) *c=0; // remove newline
      j = strlen(line);
      while (line[--j] == ' ') line[j] = 0; // remove trailing spaces
      strcat(ptitle, line+3);
      
    } else if (strstr(line, "mass")) {
      if (1 == sscanf(line, "massgoal %f", &f)) massgoal = f;
      if (1 == sscanf(line, "dmass %f", &f)) dmass = f;
    } else if (strstr(line, "fwhm")) {
      if (1 == sscanf(line, "fwhmgoal %f", &f)) fwhmgoal = f;
      if (1 == sscanf(line, "dfwhm %f", &f)) dfwhm = f;
    } else if (strstr(line, "aoe")) {
      if (1 == sscanf(line, "aoegoal %f", &f)) aoegoal = f;
      if (1 == sscanf(line, "daoe %f", &f)) daoe = f;
    }
  }
  fclose(file_in);
  printf(" Performance goals and sensitivity delats:\n"
         "  Mass: %.0f %0.f;  FWHM: %.2f %.2f; AoA)bad_fraction: %.2f %.2f\n",
         massgoal, dmass, fwhmgoal, dfwhm, aoegoal, daoe);


  //open output file
  file_out = fopen("opt.out", "w");

  if (ng==0) {
    printf(" ng = %d; No geometery parameters to scan!\n", ng);
    return -1;
  }
  printf(" >>> ng = %d\n flag   min   max  incr\n", ng);
  for (i=0; i<ng; i++) {
    printf(" -g%1d  %5.1f %5.1f %5.1f\n", g[i], flo[i], fhi[i], fincr[i]);
    fnow[i] = flo[i];
  }

  // now do scanning
  
  while (1) {
    // make system command
    strncpy(syscommand, command, 256);
    for (i=0; i<ng; i++) {
      sprintf(line, " -g%1d  %.1f", g[i], fnow[i]);
      strcat(syscommand, line);
    }
    strcat(syscommand, " > j.out");
    printf("%s\n", syscommand);
    // execute
    system(syscommand);

    // search for results
    fprintf(file_out, " <<<");
    for (i=0; i<ng; i++) fprintf(file_out, " -g%1d  %.1f", g[i], fnow[i]);
    fprintf(file_out, "\n"); fflush(file_out);

    system("grep Estimated j.out >> opt.out");
    system("grep Full j.out >> opt.out");
    system("grep overbias j.out >> opt.out");

    // calculate mass and DTC
    if ((c = strstr(syscommand, "icpc_fieldgen_opt"))) {
      sprintf(c, "mass     %s", c+18);
      system(syscommand);
      system("grep enriched j.out >> opt.out");
      // printf(" > %s\n", syscommand);
      if (do_dtc) {
        if (do_dtc > 1) {
          sprintf(c, "dtc -dd %s", c+8);
        } else {
          sprintf(c, "dtc %s", c+5);
        }
        // printf("%s\n", syscommand);
        system(syscommand);
        system("grep \" 0.1% threshold\" j.out >> opt.out");
        system("grep \" A/E fraction\" j.out >> opt.out");
        // printf(" >> %s\n", syscommand);
      }
    }

    fseek(file_out, 0, SEEK_END);
    fprintf(file_out, "\n"); fflush(file_out);

    // increment geometry values
    j = ng;
    while (j >= 0) {
      fnow[j] += fincr[j];
      if (fnow[j] < fhi[j] + fincr[j]/2.0) break;
      fnow[j] = flo[j];
      j--;
    }
    if (j < 0) {  // finished scan
      fclose(file_out);
      break;
    }
  }

  // now parse the ouput file and summarize results in a new output
  fclose (file_out);
  file_in  = fopen("opt.out", "r");
  file_out = fopen("opt2.out", "w");
  f0 = flo[ng-1];
  printf("ng = %d, g[ng-2] = %d, fo = %.1f\n", ng, g[ng-2], f0);
  fprintf(file_out, "# command: %s\n#", command);
  for (i=0; i<ng; i++) {
    fprintf(file_out, " %s ", title[g[i]]);
  }
  fprintf(file_out, " v_pinch v_depl    E_min     FOM      mass   FOM_mass   FOM_tot");
  if (do_dtc)
    fprintf(file_out, "  0.1%%_WP     FWHM_impact       A/E_fractions    FOM_final ");
  fprintf(file_out, "  E_min at (r,z)");

  float dE = 10.0, dV = 200.0, curve = 2000.0/(dE*dV);
  while (fgets(line, 256, file_in)) {
    if ((c = strstr(line, " <<< -g"))) {
      while ((c = strstr(c+2, "-g")) && 
             (sscanf(c, "-g%d %f", &j, &f) == 2)) {
        fnow[j] = f;
      }
      if (fnow[g[ng-1]] <= f0) {
        fprintf(file_out, "\n");
        f0 = fnow[g[ng-1]];
      }
      for (i=0; i<ng; i++) {
        fprintf(file_out, "%7.1f ", fnow[g[i]]);
      }
    } else if ((c = strstr(line, "Estimated depletion voltage"))) {
      c = strstr(c, "= ");
      sscanf(c+2, "%f", &f);
      if (Vhigh < f) Vhigh = f;
      if (Vlow  > f) Vlow = f;
      x = (4000.0 - f)/dV;
      fprintf(file_out, "%8.0f %6.0f", 0.0, f);
    } else if ((c = strstr(line, "Estimated pinch-off voltage"))) {
      c = strstr(c, "= ");
      sscanf(c+2, "%f", &f);
      fprintf(file_out, "%8.0f ", f);
    } else if ((c = strstr(line, "Full depletion (max pinch-off voltage)"))) {
      c = strstr(c, "= ");
      sscanf(c+2, "%f", &f);
      if (Vhigh < f) Vhigh = f;
      if (Vlow  > f) Vlow  = f;
      x = (4000.0 - f)/dV;
      fprintf(file_out, "%6.0f", f);
    } else if ((c = strstr(line, "Minimum bulk field at"))) {
      c = strstr(c, "overbias) is");
      sscanf(c+12, "%f", &f);
      if (Ehigh < f) Ehigh = f;
      if (Elow  > f) Elow  = f;
      y = (f - 200.0)/dE;
      fprintf(file_out, "%10.2f", f);

      /* calculate figure of merit
         For Emin and Vdepl:
         Emin-200 = 2000/(4000-Vdepl - dV*fom) + dE*fom
            steps of dV = 10 V/cm in Emin and dV = 200 V in Vdepl
         x = (4000 - Vdepl)/dV;   y = (Emin - 200)/dE
         a*fom^2 + b*fom + c = 0   =>   a = 1;  b = -(x'+ y');  c = x'*y' - 2000/(dE*dV)
         Then add 1 to FOM

         Incluse mass, WPfrac, FWHM, and A/E_bad_fraction in the FOM...
           FOM_mass = (mass-massgoal)/dmass
           FOM += FOM_mass
         if -d or -dd option (include DTC) then
           FOM += (wpfrac-60.0)/40.0;
           FOM += (fwhmgoal - fwhm2)/dfwhm;
           Final FOM += (aoegoal - aoe2)/daoe;

         goal snd sensitivity paramaeters can be modified in opt.config file
       */
      A = 1.0;
      B = -(x + y);
      C = (x*y - curve);
      fom = 1 + (-B - sqrt(B*B - 4.0*A*C))/(2.0*A);
      if (FOMhigh < fom) FOMhigh = fom;
      if (FOMlow  > fom) FOMlow  = fom;
      if (max_fom < fom) {
        max_fom = fom;
        for (j=0; j<8; j++) max_f[j] = fnow[j];
      }
      fprintf(file_out, "%9.3f   ", fom);
      // save (r,z)of location for E_min
      strncpy(emin_loc, strstr(c, "(r,z) = ") + 8, sizeof(emin_loc));
    } else if ((c = strstr(line, "final   mass"))) {
      c = strstr(c, "[");
      sscanf(c+1, "%f", &f);
      if (Masshigh < f) Masshigh = f;
      if (Masslow  > f) Masslow  = f;
      fomm = (f-massgoal)/dmass;
      fomt = fom + fomm;
      fprintf(file_out, "  %4.0f %9.3f %9.3f", f, fomm, fomt);
      if (!do_dtc) {
        if (FOMhigh < fomt) FOMhigh = fomt;
        if (FOMlow  > fomt) FOMlow  = fomt;
        wpfrac = fwhm1 = fwhm2 = 0;
        FWHMhigh = FWHMlow+1;
        AoEhigh = AoElow+1;
        fprintf(file_out, "      %4.1f      %.3f %.3f 0 0 0 %8.3f  %s", wpfrac, fwhm1, fwhm2, fomt, emin_loc);
      }
    } else if ((c = strstr(line, "% of non-zero WP"))) {
      c = strstr(line, " (");
      sscanf(c+2, "%f", &wpfrac);
      if (WPhigh < wpfrac) WPhigh = wpfrac;
      if (WPlow  > wpfrac) WPlow  = wpfrac;
      fprintf(file_out, "   %5.1f", wpfrac);
      fomt += (wpfrac-60.0)/40.0;
      //fprintf(file_out, " %.3f", (wpfrac-60.0)/40.0);
    } else if ((c = strstr(line, "threshold DTC"))) {
      sscanf(line, "%f %f %f", &f, &fwhm1, &fwhm2);
      fprintf(file_out, "       %.3f %.3f", fwhm1, fwhm2);
      if (FWHMhigh < fwhm2) FWHMhigh = fwhm2;
      if (FWHMlow  > fwhm2) FWHMlow  = fwhm2;
      fomt += (fwhmgoal - fwhm2)/dfwhm;
      //fprintf(file_out, " %.3f", (fwhmgoal - fwhm2)/dfwhm);
    } else if ((c = strstr(line, "A/E fraction"))) {
      sscanf(line, "%f%% %f%%", &aoe1, &aoe2);
      fprintf(file_out, "    %.3f %.3f", aoe1, aoe2);
      aoe2 += aoe1;
      fprintf(file_out, " %.3f", aoe2);
      if (AoEhigh < aoe2) AoEhigh = aoe2;
      if (AoElow  > aoe2) AoElow  = aoe2;
      if (FOMhigh < fomt) FOMhigh = fomt;
      if (FOMlow  > fomt) FOMlow  = fomt;
      fomt += (aoegoal - aoe2)/daoe;
      //fprintf(file_out, " %.3f", (aoegoal - aoe2)/daoe);
      fprintf(file_out, "%9.3f    %s", fomt, emin_loc);
    }
  }
  fprintf(file_out, "\n");
  for (i=0; i< 15; i++) fprintf(file_out, "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0\n\n");
  fclose(file_in);
  fclose(file_out);
  printf("\n Maximum FOM value is %.3f\n at", max_fom);
  for (i=0; i<ng; i++) printf(" %s ", title[g[i]]);
  printf("\n  ");
  for (i=0; i<ng; i++) printf("%7.1f ", max_f[g[i]]);
  printf("\n\n");

  /* create a stream editor that will change opt2_template.pdc to appropriate ranges */
  printf("FWHM range %.2f to %.2f\n", FWHMlow, FWHMhigh);
  f = Vhigh - Vlow; Vlow -= f/10.0; Vhigh += f/10.0;
  f = Ehigh - Elow; Elow -= f/10.0; Ehigh += f/10.0;
  f = FOMhigh - FOMlow; FOMlow -= f/10.0; FOMhigh += f/10.0;
  f = Masshigh - Masslow; Masslow -= f/10.0; Masshigh += f/10.0;
  f = WPhigh - WPlow; WPlow -= f/10.0; WPhigh += f/10.0;
  f = FWHMhigh - FWHMlow; FWHMlow -= f/10.0; FWHMhigh += f/10.0;
  f = AoEhigh - AoElow; AoElow -= f/10.0; AoEhigh += f/10.0;
  // printf("FWHM range %.2f to %.2f\n", FWHMlow, FWHMhigh);
  if (Vhigh > 3700 && Vhigh < 4010) Vhigh = 4010; 
  if (Ehigh < 220) Ehigh = 220; 
  if (Elow  > 190) Elow  = 190; 
  if (FOMlow > 0)  FOMlow  = 0; 
  if (FOMlow < -1) FOMlow  = -1;
  if (AoElow < 0)  AoElow  = 0;
 
  file_out = fopen("sed.sh", "w");
  fprintf(file_out,
          "cp opt2_template.pdc j.pdc\n"
          "sed -e 's@Config file@%s@g' -i bak j.pdc\n"
          "sed -e 's/AoElow/%.2f/g' -e 's/DeltaAoE/%.2f/g' -i bak j.pdc\n"
          "sed -e 's/Vlow/%.0f/g' -e 's/DeltaV/%.0f/g' -e 's/Vhigh/%.0f/g' -i bak j.pdc\n"
          "sed -e 's/Elow/%.0f/g' -e 's/DeltaE/%.0f/g' -e 's/Ehigh/%.0f/g' -i bak j.pdc\n"
          "sed -e 's/FOMlow/%.1f/g' -e 's/DeltaFOM/%.1f/g' -e 's/FOMhigh/%.1f/g' -i bak j.pdc\n"
          "sed -e 's/Masslow/%.0f/g' -e 's/DeltaMass/%.0f/g' -e 's/WPlow/%.0f/g' -e 's/DeltaWP/%.0f/g' -i bak j.pdc\n"
          "sed -e 's/FWHMlow/%.2f/g' -e 's/DeltaFWHM/%.2f/g' -i bak j.pdc\n"
          "mv j.pdc opt2.pdc\n",
          ptitle, AoElow, AoEhigh-AoElow, Vlow, Vhigh-Vlow, Vhigh-1, Elow, Ehigh-Elow, Ehigh-1, FOMlow,
          FOMhigh-FOMlow, FOMhigh-0.1, Masslow, Masshigh-Masslow, WPlow, WPhigh-WPlow, FWHMlow, FWHMhigh-FWHMlow);
  fclose(file_out);
  system("source sed.sh");
    //system("source sed.sh; rm sed.sh");

  return 0;
}


