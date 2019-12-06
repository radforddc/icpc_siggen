/* NB: This code is not as flexible as it could be. Non-integer dimensions
   may not be handled properly ("implemented" but not tested)
   Also: Add setup_done flag?
   Modified November 2007: Asymmetric hexagon, two different detector shapes/KL
   Also: setup mode, to compensate for deficiensies in field tables
   setup mode affects in_crystal and its helper function project_to_edges
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "detector_geometry.h"
#include "point.h"
#include "calc_signal.h"
#include "signal_calc_util.h"

float zmax_detector(MJD_Siggen_Setup *setup){
  return setup->xtal_length;
}

float rmax_detector(MJD_Siggen_Setup *setup){
  return setup->xtal_radius;
}

float hole_r_detector(MJD_Siggen_Setup *setup){
  return setup->hole_radius;
}

float hole_z_detector(MJD_Siggen_Setup *setup){
  return setup->hole_gap;
}


#define SQ(x) ((x)*(x))

/* returns 0 (false) or 1 (true) depending on whether pt is inside the crystal */
int in_crystal(MJD_Siggen_Setup *setup, point pt) {
  float z, r, rsq, r2, rmax, zmin, zmax, dx, dy;
  static int first = 1, hole_offset = 0, pc_offset = 0;
  float lith  = setup->Li_thickness;
  float r1, z1, br, a, b;

  if (first) {
    first = 0;
    setup->pc_offset_max = 0;
    if (setup->pc_offset_x !=0 ||
        setup->pc_offset_y !=0 ) {
      pc_offset = 1;
      setup->pc_offset_max = sqrt(SQ(setup->pc_offset_x) + SQ(setup->pc_offset_y));
    }
    setup->hole_offset_max = 0;
    if (setup->hole_offset_x_top !=0 ||
        setup->hole_offset_y_top !=0 ||
        setup->hole_offset_x_bottom !=0 ||
        setup->hole_offset_y_bottom !=0) {
      hole_offset = 1;
      setup->hole_offset_max = sqrt(SQ(setup->hole_offset_x_top) + SQ(setup->hole_offset_y_top));
      r2 = sqrt(SQ(setup->hole_offset_x_bottom) + SQ(setup->hole_offset_y_bottom));
      if (setup->hole_offset_max < r2) setup->hole_offset_max = r2;
    }
  }

  zmin = 0;
  zmax = setup->xtal_length - lith;
  rmax = setup->xtal_radius - lith;

  z = pt.z;
  if (z >= zmax || z < 0) return 0;
  rsq = SQ(pt.x) + SQ(pt.y);
  if (rsq > SQ(rmax)) return 0;

  r = sqrt(rsq);
  r1 = rmax - r;  // distance from outer radius
  z1 = zmax - z;  // distance from top of crystal

  if (pc_offset) {
    r2 = sqrt(SQ(pt.x - setup->pc_offset_x) + SQ(pt.y - setup->pc_offset_y));
  } else {
    r2 = r;
  }
  if (setup->ditch_depth > 0 && r2 >= setup->wrap_around_radius && z < lith) return 0;
  /* check point contact */
  if (z <= setup->pc_length && r2 <= setup->pc_radius) {
    if (!setup->bulletize_PC) return 0;
    if (setup->pc_length > setup->pc_radius) {
      a = setup->pc_length - setup->pc_radius;
      if (z <= a || SQ(z-a) + SQ(r2) <= SQ(setup->pc_radius)) return 0;
    } else {
      a = setup->pc_radius - setup->pc_length;
      if (r2 <= a || SQ(z) + SQ(r2-a) <= SQ(setup->pc_length)) return 0;
    }
  }
  /* check ditch */
  if (z <= setup->ditch_depth  &&
      setup->ditch_thickness > 0 && setup->wrap_around_radius > 0 &&
      r2 <= setup->wrap_around_radius &&
      r2 >= setup->wrap_around_radius - setup->ditch_thickness) return 0;

  /* check central hole */
  if (hole_offset) {
    dx = setup->hole_offset_x_top +
      (setup->hole_offset_x_bottom - setup->hole_offset_x_top) *
      (setup->xtal_length - pt.z) / setup->hole_length;
    dy = setup->hole_offset_y_top +
      (setup->hole_offset_y_bottom - setup->hole_offset_y_top) *
      (setup->xtal_length - pt.z) / setup->hole_length;
    r2 = sqrt(SQ(pt.x - dx) + SQ(pt.y - dy));
  } else {
    r2 = r;
  }
  if (r2 <= setup->hole_radius + lith &&
      z1 <= setup->hole_length) {  // note no lith added here since it was already subtracted from zmax
    br = setup->hole_bullet_radius + 0.8*lith;
    b = zmax - setup->hole_length + br;
    if (z >= b) return 0;
    a = setup->hole_radius + lith - br;
    if (r2 <= a || SQ(b-z) + SQ(a-r2) <= SQ(br)) return 0;
  }
  /* check inner taper of hole */
  if (z1 <= setup->inner_taper_length &&
      r2 <= setup->hole_radius + lith +
            ((setup->inner_taper_length - z1) *
              setup->inner_taper_width / setup->inner_taper_length)) return 0;      

  /* check outer taper of crystal */
  if (z1 <= setup->outer_taper_length &&
      r1 <= ((setup->outer_taper_length - z1) *
              setup->outer_taper_width / setup->outer_taper_length)) return 0;
  /* check 45-degree bottom outer taper of crystal */
  if ( z <= setup->bottom_taper_length + 0.71*lith &&
      r1 <= setup->bottom_taper_length + 0.71*lith - z) return 0;

  /* check bulletizations */
  br = setup->top_bullet_radius - lith * 0.7;
  if (z1 <= br &&
      r1 <= br &&
      SQ(br - r1) + SQ(br - z1) >= br*br) return 0;
  br = setup->bottom_bullet_radius - lith * 0.7;
  if ( z <= br + lith &&
      r1 <= br &&
      SQ(br - r1) + SQ(br - z + lith) >= br*br) return 0;

  return 1;
}
#undef SQ
