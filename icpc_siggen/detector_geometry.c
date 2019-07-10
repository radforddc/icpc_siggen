/* detector_geometry_ppc.c -- for "ppc" geometry
 * Karin Lagergren
 *
 * This module keeps track of the detector geometry
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "detector_geometry.h"
#include "point.h"
#include "cyl_point.h"

#define SQ(x) ((x)*(x))
/* outside_detector
   returns 1 if pt is outside the detector, 0 if inside detector
*/
int outside_detector(point pt, MJD_Siggen_Setup *setup){
  float r, z, br, a, b;

  z = pt.z;
  if (z > setup->zmax || z < 0) return 1;

  r = sqrt(SQ(pt.x)+SQ(pt.y));
  if (r > setup->rmax) return 1;
  br = setup->top_bullet_radius;
  if (z > setup->zmax - br &&
      r > (setup->rmax - br) + sqrt(SQ(br)- SQ(z-(setup->zmax - br)))) return 1;
  if (setup->pc_radius > 0 &&
      z < setup->pc_length && r < setup->pc_radius) {
    if (!setup->bulletize_PC) return 1;
    if (setup->pc_length > setup->pc_radius) {
      a = setup->pc_length - setup->pc_radius;
      if (z < a || SQ(z-a) + SQ(r) < SQ(setup->pc_radius)) return 1;
    } else {
      a = setup->pc_radius - setup->pc_length;
      if (r < a || SQ(z) + SQ(r-a) < SQ(setup->pc_length)) return 1;
    }
    return 0;
  }
  if (setup->bottom_taper_length > 0 && z < setup->bottom_taper_length &&
      r > setup->zmax - setup->bottom_taper_length + z) return 1;
  if (setup->ditch_depth > 0 && z < setup->ditch_depth  &&
      setup->ditch_thickness > 0 && setup->wrap_around_radius > 0 &&
      r < setup->wrap_around_radius &&
      r > setup->wrap_around_radius - setup->ditch_thickness) return 1;

  /* check hole */
  if (r < setup->hole_radius &&
      z > setup->zmax - setup->hole_length) {
    b = setup->zmax - setup->hole_length + setup->hole_bullet_radius;
    if (z > b) return 1;
    a = setup->hole_radius - setup->hole_bullet_radius;
    if (r < a || SQ(b-z) + SQ(r-a) < SQ(setup->hole_bullet_radius)) return 1;
    return 0;
  }

  /* check outer taper of crystal */
  if (setup->outer_taper_length > 0 &&
      z > setup->zmax - setup->outer_taper_length &&
      r > setup->rmax - ((z - setup->zmax + setup->outer_taper_length) *
                         setup->outer_taper_width / setup->outer_taper_length)) return 1;
  /* check inner taper of hole */
  if (setup->inner_taper_length > 0 &&
      z > setup->zmax - setup->inner_taper_length &&
      r < setup->hole_radius + ((z - setup->zmax + setup->inner_taper_length) *
                         setup->inner_taper_width / setup->inner_taper_length)) return 1;

  return 0;
}

int outside_detector_cyl(cyl_pt pt, MJD_Siggen_Setup *setup){
  float r, z, br, a, b;

  z = pt.z;
  if (z > setup->zmax || z < 0) return 1;

  r = pt.r;
  if (r > setup->rmax) return 1;
  br = setup->top_bullet_radius;
  if (z > setup->zmax - br &&
      r > (setup->rmax - br) + sqrt(SQ(br)- SQ(z-(setup->zmax - br)))) return 1;
  if (setup->pc_radius > 0 &&
      z < setup->pc_length && r < setup->pc_radius) {
    if (!setup->bulletize_PC) return 1;
    if (setup->pc_length > setup->pc_radius) {
      a = setup->pc_length - setup->pc_radius;
      if (z < a || SQ(z-a) + SQ(r) < SQ(setup->pc_radius)) return 1;
    } else {
      a = setup->pc_radius - setup->pc_length;
      if (r < a || SQ(z) + SQ(r-a) < SQ(setup->pc_length)) return 1;
    }
    return 0;
  }
  if (setup->bottom_taper_length > 0 && z < setup->bottom_taper_length &&
      r > setup->zmax - setup->bottom_taper_length + z) return 1;
  if (setup->ditch_depth > 0 && z < setup->ditch_depth  &&
      setup->ditch_thickness > 0 && setup->wrap_around_radius > 0 &&
      r < setup->wrap_around_radius &&
      r > setup->wrap_around_radius - setup->ditch_thickness) return 1;

  /* check hole */
  if (r < setup->hole_radius &&
      z > setup->zmax - setup->hole_length) {
    b = setup->zmax - setup->hole_length + setup->hole_bullet_radius;
    if (z > b) return 1;
    a = setup->hole_radius - setup->hole_bullet_radius;
    if (r < a || SQ(b-z) + SQ(r-a) < SQ(setup->hole_bullet_radius)) return 1;
  }
      
  /* check outer taper of crystal */
  if (setup->outer_taper_length > 0 &&
      z > setup->zmax - setup->outer_taper_length &&
      r > setup->rmax - ((z - setup->zmax + setup->outer_taper_length) *
                         setup->outer_taper_width / setup->outer_taper_length)) return 1;
  /* check inner taper of hole */
  if (setup->inner_taper_length > 0 &&
      z > setup->zmax - setup->inner_taper_length &&
      r < setup->hole_radius + ((z - setup->zmax + setup->inner_taper_length) *
                         setup->inner_taper_width / setup->inner_taper_length)) return 1;

  return 0;
}
#undef SQ
