/* cyl_point.c
 * Karin Lagergren
 * simple routines for handling cylindrical coordinates
 */
#include <stdio.h>
#include <math.h>

#include "cyl_point.h"
#include "point.h"


struct cyl_pt cart_to_cyl(struct point cart_pt){
  struct cyl_pt cyl_pt;

  cyl_pt.r = sqrt(cart_pt.x*cart_pt.x + cart_pt.y*cart_pt.y);
  if (cart_pt.x == 0.0){
    if (cart_pt.y > 0.0) cyl_pt.phi = M_PI/2;
    else cyl_pt.phi = -M_PI/2;
  }else{
    cyl_pt.phi = atan(cart_pt.y/cart_pt.x);
    if (cart_pt.x < 0) cyl_pt.phi += M_PI;
  }
  cyl_pt.z = cart_pt.z;
  
  return cyl_pt;
}

struct point cyl_to_cart(struct cyl_pt cyl_pt){
  struct point cart_pt;

  cart_pt.x = cyl_pt.r * cos(cyl_pt.phi);
  cart_pt.y = cyl_pt.r * sin(cyl_pt.phi);
  cart_pt.z = cyl_pt.z;

  return cart_pt;
}

#define SQUARED(x) ((x)*(x))
float cart_distance(struct point pt1, struct point pt2){
  return sqrt(SQUARED(pt1.x-pt2.x) + SQUARED(pt1.y-pt2.y) 
	      + SQUARED(pt1.z-pt2.z));
}
#undef SQUARED


float cyl_distance(struct cyl_pt pt1, struct cyl_pt pt2){
  /*yes, I'm very lazy. It's one of my best traits*/
  return cart_distance(cyl_to_cart(pt1), cyl_to_cart(pt2));
}

struct cyl_pt cyl_diff(struct cyl_pt from, struct cyl_pt to){
  struct cyl_pt v;

  v.r = from.r - to.r;
  v.phi = from.phi - to.phi;
  v.z = from.z - to.z;

  return v;
}

float vector_norm_cyl(struct cyl_pt pt, struct cyl_pt *e){
  float norm;

  norm = sqrt(pt.r*pt.r + pt.z*pt.z);
  e->phi = pt.phi;
  e->r = pt.r/norm;
  e->z = pt.z/norm;

  return norm;
}


