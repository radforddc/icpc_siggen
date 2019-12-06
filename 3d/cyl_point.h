/* cyl_point.h
 * Karin Lagergren
 * simple routines for handling cylindrical coordinates
 */
#ifndef _CYL_POINT_H
#define _CYL_POINT_H

#include "point.h"

struct cyl_pt{
  float r;
  float phi;
  float z;
};

typedef struct cyl_pt cyl_pt;

struct cyl_int_pt{
  int r;
  int phi;
  int z;
};

typedef struct cyl_int_pt cyl_int_pt;

struct cyl_pt cart_to_cyl(struct point cart_pt);
struct point cyl_to_cart(struct cyl_pt cyl_pt);

float cart_distance(struct point pt1, struct point pt2);
float cyl_distance(struct cyl_pt pt1, struct cyl_pt pt2);

struct cyl_pt cyl_diff(struct cyl_pt from, struct cyl_pt to);

float vector_norm_cyl(struct cyl_pt pt, struct cyl_pt *e);


#endif /*#ifndef _CYL_POINT_H*/
