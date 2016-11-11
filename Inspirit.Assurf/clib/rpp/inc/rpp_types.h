#ifndef __RPP_TYPES_H__
#define __RPP_TYPES_H__

#include "rpp_const.h"

typedef int bool;
typedef double rpp_float;
typedef double real_t;
typedef double rpp_vec[3];
typedef double rpp_mat[9];

// standard types
//
typedef real_t vec3[3];
typedef struct { vec3 v; }vec3_t;
typedef real_t mat33[9];
typedef struct { mat33 m; }mat33_t;

typedef struct
{
	mat33_t R;
	vec3_t  t;
	real_t  E;
	mat33_t PoseLu_R;
	vec3_t  PoseLu_t;
	real_t  obj_err;
}pose_t;

typedef struct
{
	vec3_t v;
	real_t  s;
}quat_t;

typedef struct
{
	mat33_t initR;
	real_t tol;
	real_t epsilon;
	unsigned int max_iter;
}options_t;

#endif
