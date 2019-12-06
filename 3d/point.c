/* point.c
 * Karin Lagergren
 * 
 * Cartesian points and vectors for signal calculation code
 *
 * November 2007, added dot product, vector subtraction, cross product 
 */

#include <stdio.h>
#include <math.h>
#include "point.h"


#define SQ(x) ((x)*(x))

/* norm
 * returns the absolute value of vector v
 * e is unit vector in v's direction
 */
float vector_norm(vector v, vector *e){
  float length;
  
  length = vector_length(v);

  e->x = v.x / length;
  e->y = v.y / length;
  e->z = v.z / length;
  
  return length;
}

/* vector_length
 * returns the length of vector v
 */
float vector_length(vector v){
  return sqrt(SQ(v.x) + SQ(v.y) + SQ(v.z));
}

/* distance
 * returns distance between points
 */
float distance(point pt1, point pt2){
  float d;

  d = sqrt(SQ(pt1.x - pt2.x) + SQ(pt1.y - pt2.y) + SQ(pt1.z - pt2.z));

  return d;
}

/* vector_add */
vector vector_add(vector v1, vector v2){
  vector v;
  
  v.x = v1.x + v2.x;
  v.y = v1.y + v2.y;
  v.z = v1.z + v2.z;

  return v;
}

/*vector subtraction, v1-v2*/
vector vector_sub(vector v1, vector v2){
  vector v;
  
  v.x = v1.x - v2.x;
  v.y = v1.y - v2.y;
  v.z = v1.z - v2.z;
  
  return v;
}

/*vector dot product, v1*v1*/
float dot_prod(vector v1, vector v2){
  return  v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
}

/*vector cross product, v1 X v2*/
vector cross_prod(vector v1, vector v2){
  vector v;
  
  v.x = v1.y*v2.z - v1.z*v2.y;
  v.y = v1.z*v2.x - v1.x*v2.z;
  v.z = v1.x*v2.y - v1.y*v2.x;

  return v;
}



/*scale vector by factor -- better name?*/
vector vector_scale(vector v, float factor){
  vector res;

  res.x = v.x*factor;
  res.y = v.y*factor;
  res.z = v.z*factor;

  return res;
}

/*rotate vector by angle_deg degrees*/
vector vector_rotate_z(vector v, float angle_deg){
  float angle_rad;
  vector res;
  angle_rad = angle_deg/180.0*M_PI;

  res.x = cos(angle_rad)*v.x - sin(angle_rad)*v.y;
  res.y = sin(angle_rad)*v.x + cos(angle_rad)*v.y;
  res.z = v.z;

  return res;
}


/* pt_to_str
 * string representation of point pt
 * written to string str, which has guaranteed length len 
 */
int pt_to_str(char *str, int len, point pt){
  snprintf(str, len, "(%.1f %.1f %.1f)", pt.x, pt.y, pt.z);
  return 0;
}

#undef SQ
