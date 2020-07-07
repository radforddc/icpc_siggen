#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include "mjd_siggen.h"

/* calculates approximate mass of detector, using 5.323 g/cm3
   does NOT take bulletizations into account!
   DCR Nov 2017

   89% enriched 76Ge density = 5.544?
   - GERDA measured 5.55  [DOI 10.1140/epjc/s10052-014-3253-0 (p3)]
 */

int main(int argc, char **argv)
{

  MJD_Siggen_Setup setup;

  /* --- default values, normally over-ridden by values in a *.conf file --- */
  /* all dimensions are in mm */
  float   r  = 0;  // radius of detector
  float   l  = 0;  // length of detector
  float   rc = 0;  // radius of central contact
  float   lc = 0;  // length of central contact
  float   lt = 0;  // length of taper
  float   ro = 0;  // radius of wrap-around outer (li) contact
  float   lo = 0;  // length of ditch next to wrap-around outer (li) contact
  float   wo = 0;  // width of ditch next to wrap-around outer (li) contact

  float   rh  = 10; // radius of core hole in outer (li) contact
  float   lh  = 80; // length of core hole in outer (li) contact
  float   hb  = 5;  // bulletization radius at bottom of hole
  float   otl = 80; // length of outer radial taper of crystal
  float   otw = 10; // width/amount of outer radial taper (decrease in radius)
  float   htl = 0;  // length of radial tapered part of core hole
  float   htw = 0;  // width/amount of radial taper (increase in radius) of hole
  float   tbr = 0;  // top bullet radius
  float   bbr = 0;  // bottom bullet radius
  /* ---  ---  ---  ---  ---  ---  ---  ---  ---  ---  ---  ---  ---  --- */

  float  v, r1, r2, h1, h2, v1, v2;
  float  pi = 3.141593;


  if (argc < 2 || read_config(argv[1], &setup)) {
    printf("Usage: %s <config_file_name>\n\n", argv[0]);
    return 0;
  }

  l  =  setup.xtal_length;
  r  =  setup.xtal_radius;
  lc =  setup.pc_length;
  rc =  setup.pc_radius;
  ro =  setup.wrap_around_radius;
  lo =  setup.ditch_depth;
  wo =  setup.ditch_thickness;
  lt =  setup.bottom_taper_length;
  lh =  setup.hole_length;
  rh =  setup.hole_radius;
  hb =  setup.hole_bullet_radius;
  otl = setup.outer_taper_length;
  otw = setup.outer_taper_width;
  htl = setup.inner_taper_length;
  htw = setup.inner_taper_width;
  tbr = setup.top_bullet_radius;
  bbr = setup.bottom_bullet_radius;


  printf("\n\n"
         "      Crystal: Radius x Length: %.1f x %.1f mm\n"
         " Bottom taper: %.1f mm\n", r, l, lt);
  printf("    Core hole: Radius x length: %.1f x %.1f mm, taper %.1f x %.1f mm (%2.f degrees),"
         " bullet radius %.1f mm\n", rh, lh, htw, htl, setup.taper_angle, hb);
  printf("Outside taper: %.1f mm x %.1f mm (%.2f degrees)\n", otw, otl, setup.taper_angle);
  printf("      Contact: Radius x length: %.1f x %.1f mm\n", rc, lc);
  printf("  Wrap-around: Radius x ditch x gap:  %.1f x %.1f x %.1f mm\n\n", ro, lo, wo);

  v  = pi * (r*r * l - rc*rc * lc);
  printf("         cylinder: %7.0f mm3  =  %4.0f g\n", v, v*0.005323);

  // hole (with optional taper)
  if (htw > 0.1) {
    r1 = rh + htw;
    r2 = rh;
    h2 = htl * rh / htw;
    h1 = h2 + htl;
    v1 = pi * r1*r1 * h1 / 3.0;
    v2 = pi * r2*r2 * h2 / 3.0;
    v -= v1 - v2;
    printf("     tapered hole: %7.0f mm3  =  %4.0f g\n", v1-v2, (v1-v2)*0.005323);
  }
  if (rh > 0 && lh > htl) {
    v1 = pi * rh*rh * (lh - htl);
    v -= v1;
    printf("    straight hole: %7.0f mm3  =  %4.0f g\n", v1, v1*0.005323);
  }

  // bottom taper
  if (lt > 0.1) {
    r1 = h1 = r;
    r2 = h2 = r - lt;
    v1 = pi * r1*r1 * h1 / 3.0;
    v2 = pi * r2*r2 * h2 / 3.0;
    v1 = pi * r*r * lt - v1 + v2;
    v -= v1;
    printf("     bottom taper: %7.0f mm3  =  %4.0f g\n", v1, v1*0.005323);
  }

  // outer (top) taper
  if (otl * otw > 0.1) {
    r1 = r;
    r2 = r - otw;
    h1 = otl * r / otw;
    h2 = h1 - otl;
    v1 = pi * r1*r1 * h1 / 3.0;
    v2 = pi * r2*r2 * h2 / 3.0;
    v1 = pi * r*r * otl - v1 + v2;
    v -= v1;
    printf("outer (top) taper: %7.0f mm3  =  %4.0f g\n", v1, v1*0.005323);
  }

  // top and bottom bullet radius
  if (tbr > 0.1) {
    v1 = 1.35 * r * tbr*tbr;
    v -= v1;
    printf("top bullet radius: %7.0f mm3  =  %4.0f g\n", v1, v1*0.005323);
  }
  if (bbr > 0.1) {
    v1 = 1.35 * r * bbr*bbr;
    v -= v1;
    printf("  bottom bullet r: %7.0f mm3  =  %4.0f g\n", v1, v1*0.005323);
  }

  // outer (top) taper
  if (otl * otw > 0.1) {
    r1 = r;
    r2 = r - otw;
    h1 = otl * r / otw;
    h2 = h1 - otl;
    v1 = pi * r1*r1 * h1 / 3.0;
    v2 = pi * r2*r2 * h2 / 3.0;
    v1 = pi * r*r * otl - v1 + v2;
    v -= v1;
    printf("outer (top) taper: %7.0f mm3  =  %4.0f g\n", v1, v1*0.005323);
  }

  // ditch
  if (ro * lo > 0.1) {
    r1 = ro;
    r2 = r1 - wo;
    v1 = pi * (r1*r1 - r2+r2) * lo;
    v -= v1;
    printf("            ditch: %7.0f mm3  =  %4.0f g\n", v1, v1*0.005323);
  }

  printf("\n  final volume:  %.0f mm3  =  %.1f cm3\n", v, v/1000.0);
  printf("  final   mass:  %.0f g   [%.0f g for enriched 76Ge]\n\n", v * 5.323/1000.0, v * 5.544/1000.0);

  /* summarize geometry and mass for enriched detectors */
  if (1) {
    char *c, name[256];
    strncpy(name, argv[1], 256);
    while ((c = strstr(name, "/"))) strncpy(name, c+1, 256);
    if ((c = strstr(name, ".con"))) *c = 0;
    printf(" detector  len  diam  enr_mass Vop  well_diamxlen taper ditchOD\n"
           " %8s %5.1f %5.1f %6.0f  %4.0f  %4.1fx%.1f   %4.1f %4.1f\n",
           name, l, 2.0*r, v * 5.544/1000.0, setup.xtal_HV, 2.0*rh, lh, htl, 2.0*ro);
  }

  return 0;
}
