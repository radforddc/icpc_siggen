#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

float fminf(float x, float y);

#include "mjd_siggen.h"

static int   in_crystal(MJD_Siggen_Setup *setup, point pt);
static int   init_point_fractions(MJD_Siggen_Setup *setup);
static float get_fraction(MJD_Siggen_Setup *setup, int i, int j, int k, int nsteps);


int geometry_init(MJD_Siggen_Setup *setup) {

  init_point_fractions(setup);
  return 1;
}

int init_ev_calc(MJD_Siggen_Setup *setup) {
  int i, j, k;
  point pt;
  float r;

  for (i = 0; i < setup->numx; i++) {
    pt.x = setup->xmin + i*setup->xtal_grid;
    for (j = 0; j < setup->numy; j++) {
      pt.y = setup->ymin + j*setup->xtal_grid;
      r = sqrt(pt.x*pt.x + pt.y*pt.y);
      for (k = 0; k < setup->numz; k++) {
	pt.z = setup->zmin + k*setup->xtal_grid;
	if (setup->point_type[i][j][k] == CONTACT_0) {
	  setup->v[0][i][j][k] = setup->v[1][i][j][k] = 0.0;
	} else if (setup->point_type[i][j][k] == CONTACT_VB) {
	  setup->v[0][i][j][k] = setup->v[1][i][j][k] = setup->xtal_HV;
	} else if (setup->point_type[i][j][k]  == OUTSIDE) {
	  setup->v[0][i][j][k] = setup->v[1][i][j][k] = 0.0;
	} else {
	  setup->v[0][i][j][k] = setup->v[1][i][j][k] = setup->xtal_HV*(setup->xtal_radius - r) /
            (setup->xtal_radius - setup->hole_radius);
	  if (setup->v[0][i][j][k] >= setup->xtal_HV*0.95) setup->v[0][i][j][k] = setup->xtal_HV*0.95;
	  if (setup->v[0][i][j][k] <= setup->xtal_HV*0.05) setup->v[0][i][j][k] = setup->xtal_HV*0.05;
	}
      }
    }
  }
  printf("ev ... init done\n"); fflush(stdout);
  return 0;

}

int init_wp_calc(MJD_Siggen_Setup *setup) {
  int i, j, k;
  point pt;
  float r;

  for (i = 0; i < setup->numx; i++) {
    for (j = 0; j < setup->numy; j++) {
      for (k = 0; k < setup->numz; k++) {
	setup->v[0][i][j][k] = setup->v[1][i][j][k] = 0.0;
      }
    }
  }

  for (i = 0; i < setup->numx; i++) {
    pt.x = setup->xmin + i*setup->xtal_grid;
    for (j = 0; j < setup->numy; j++) {
      pt.y = setup->ymin + j*setup->xtal_grid;
      r = sqrt(pt.x*pt.x + pt.y*pt.y);
      for (k = 0; k < setup->numz; k++) {
	pt.z = setup->zmin + k*setup->xtal_grid;
	if (setup->point_type[i][j][k] == CONTACT_VB) {
          setup->v[0][i][j][k] = setup->v[1][i][j][k] = 0.0;
	} else if (setup->point_type[i][j][k] == CONTACT_0) {
          setup->v[0][i][j][k] = setup->v[1][i][j][k] = 1.0;
	} else if (setup->point_type[i][j][k]  == OUTSIDE) {
	  setup->v[0][i][j][k] = setup->v[1][i][j][k] = 0.0;
	} else {
	  setup->v[0][i][j][k] = setup->v[1][i][j][k] = 0.5;
	}
      }
    }
  }
  printf("wp ... init done\n"); fflush(stdout);
  return 0;


}

#define SQ(x) ((x)*(x))

static int init_point_fractions(MJD_Siggen_Setup *setup) {
  int   i, j, k;
  point pt;
  float z, r, f;

  printf("\n  Initializing points in crystal for x = %.1f to %.1f, grid %.2f...\n",
         setup->xmin, setup->xmax, setup->xtal_grid);
  if (setup->ditch_thickness <=0 &&
      setup->bottom_taper_length > 0) {
    if (setup->wrap_around_radius == 0 ||
        setup->wrap_around_radius > setup->xtal_radius -
        setup->bottom_taper_length - 1.71*setup->Li_thickness)
      setup->wrap_around_radius = setup->xtal_radius -
        setup->bottom_taper_length - 1.71*setup->Li_thickness;

  }

  for (i = 0; i < setup->numx; i++) {
    printf("\r %d/%d", i, setup->numx-1);fflush(stdout);
    pt.x = setup->xmin + i*setup->xtal_grid;
    for (j = 0; j < setup->numy; j++) {
      pt.y = setup->ymin + j*setup->xtal_grid;
      r = sqrt(SQ(pt.x - setup->pc_offset_x) + SQ(pt.y - setup->pc_offset_y));
      for (k = 0; k < setup->numz; k++) {
	z = pt.z = setup->zmin + k*setup->xtal_grid;
	f = get_fraction(setup, i, j, k, 3);
	if (f == 1.0) {                           // voxel is all inside bulk
	  setup->point_type[i][j][k] = INSIDE;
	} else if (f == 0.0) {                    // voxel is all outside
	  if (setup->ditch_depth > 0 &&
              z <= setup->ditch_depth &&
              r <  setup->wrap_around_radius) {                                  /* ditch */
            setup->point_type[i][j][k] = DITCH;
            setup->epsilon[i][j][k] = 1.0;
          } else {
            setup->point_type[i][j][k] = OUTSIDE;
          }
	} else {                                  // voxel is on a boundary
	  if (z < setup->pc_length + setup->xtal_grid &&
              r < setup->pc_radius + setup->xtal_grid) {                        /* point contact */
	    setup->point_type[i][j][k] = CONTACT_0;
	  } else if (k == 0 &&
                     r < setup->wrap_around_radius - setup->ditch_thickness) {  /* passivated surface */
            setup->point_type[i][j][k] = PASSIVE;
	  } else if (setup->ditch_depth > 0 &&
                     z <= setup->ditch_depth + setup->xtal_grid &&
                     r <= setup->wrap_around_radius) {                          /* ditch */
            setup->point_type[i][j][k] = DITCH;
	  } else {                                                              /* outside contact */
	    setup->point_type[i][j][k] = CONTACT_VB;
	  }
	}
      }
    }
  }

  printf("\n");
  return 0;
}

static float get_fraction(MJD_Siggen_Setup *setup, int i, int j, int k, int nsteps) {
  int ii, jj, kk;
  point pt;
  int n, m;

  // first try all eight corners of voxel to see if they are inside the detector
  n = m = 0;
  for (ii = 0; ii < 2; ii++) {
    pt.x = setup->xmin + (i + ii - 0.5)*setup->xtal_grid;
    if (pt.x < setup->xmin || pt.x > setup->xmax) {
      n += 4;
      continue;
    }
    for (jj = 0; jj < 2; jj++) {
      pt.y = setup->ymin + (j + jj - 0.5)*setup->xtal_grid;
      if (pt.y < setup->ymin || pt.y > setup->ymax) {
	n += 2;
	continue;
      }
      for (kk = 0; kk < 2; kk++) {
        n++;
	pt.z = setup->zmin + (k + kk - 0.5)*setup->xtal_grid;
	if (pt.z < setup->zmin || pt.z > setup->zmax) continue;
	if (in_crystal(setup, pt)) m++;
      }
    }
  }
  if (m == 0) return 0.0;  // all are outside
  if (m == n) return 1.0;  // all are inside

  /*
  // no definitive answer, so use (nsteps+1)^3 smaller fractions of the voxel
  n = m = 0;
  for (ii = -nsteps; ii <= nsteps; ii++) {
    pt.x = setup->xmin + (i + ii/2.0/nsteps)*setup->xtal_grid;
    if (pt.x < setup->xmin || pt.x > setup->xmax) {
      n += (2*nsteps+1) * (2*nsteps+1);
      continue;
    }
    for (jj = -nsteps; jj <=nsteps; jj++) {
      pt.y = setup->ymin + (j + jj/2.0/nsteps)*setup->xtal_grid;
      if (pt.y < setup->ymin || pt.y > setup->ymax) {
	n += 2*nsteps+1;
	continue;
      }
      for (kk = -nsteps; kk <= nsteps; kk++) {
        n++;
	pt.z = setup->zmin + (k + kk/2.0/nsteps)*setup->xtal_grid;
	if (pt.z < setup->zmin || pt.z > setup->zmax) continue;
	if (in_crystal(setup, pt)) m++;
      }
    }
  }
  */
  return ((float) m) / ((float) n);
}


/* returns 0 (false) or 1 (true) depending on whether pt is inside the crystal */
static int in_crystal(MJD_Siggen_Setup *setup, point pt) {
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

  // special hack for 9-mm notch on side of Mirion detector 60A
  // if (pt.y > setup->xtal_radius + pt.z - 9.0) return 0;
  // special hack for flat on side of ORTEC PPC detector
  // if (pt.y > setup->xtal_radius - 10.0) return 0;

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
      if (r <= a || SQ(z) + SQ(r2-a) <= SQ(setup->pc_length)) return 0;
    }
    return 1;
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

