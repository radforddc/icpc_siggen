/* point.h
 * Karin Lagergren
 * 
 * Cartesian points and vectors for signal calculation code
 *
 * November 2007, added dot product, vector subtraction, cross product
 */
#ifndef _POINT_H
#define _POINT_H

struct point{
  float x;
  float y;
  float z;
};

typedef struct point point;
typedef struct point vector;

typedef struct{
  int x;
  int y;
  int z;
} int_pt;

/* vector_norm
 * returns the absolute value of vector v
 * e is unit vector in v's direction
 */
float vector_norm(vector v, vector *e);


/* vector_length
 * returns length of vector v
 */
float vector_length(vector v);


/* distance
 * returns distance between points
 */
float distance(point pt1, point pt2);

/* vector_add */
vector vector_add(vector v1, vector v2);

/*vector subtraction, v1-v2*/
vector vector_sub(vector v1, vector v2);

/*vector dot product, v1*v1*/
float dot_prod(vector v1, vector v2);

/*vector cross product, v1 X v2*/
vector cross_prod(vector v1, vector v2);

/*scale vector by factor*/
vector vector_scale(vector v, float factor);

/*rotate vector by angle_deg degrees around z axis*/
vector vector_rotate_z(vector v, float angle_deg);

/* pt_to_str
 * string representation of point pt
 * written to string str, which has guaranteed length len 
 */
int pt_to_str(char *str, int len, point pt);


#endif /*#ifndef _POINT_H*/
