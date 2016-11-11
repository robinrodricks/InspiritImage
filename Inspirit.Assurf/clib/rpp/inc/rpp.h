#ifndef __RPP_H__
#define __RPP_H__

#include "rpp_types.h"

//namespace rpp {
// ------------------------------------------------------------------------------------------
void Quaternion_byAngleAndVector(quat_t *Q, const real_t q_angle, const vec3_t *q_vector);
void GetRotationbyVector(mat33_t *R, const vec3_t *v1, const vec3_t *v2);
void xform(vec3_t* Q, const vec3_t *P, const mat33_t *R, const vec3_t *t, const int n);
void xformproj(vec3_t *Qp, const vec3_t *P, const mat33_t *R, const vec3_t *t, const int n);
void rpyMat(mat33_t *R, const vec3_t *rpy);
void rpyAng(vec3_t *angs, const mat33_t *R);
void rpyAng_X(vec3_t *ang_zyx, const mat33_t *R);
void decomposeR(mat33_t *Rz, const mat33_t *R);
void optimal_t(vec3_t *t ,const mat33_t *R, const mat33_t *G, const mat33_t *F,const vec3_t *P, int n);
void abskernel(mat33_t *R, vec3_t *t, vec3_t *Qout, real_t *err2, 
			   const vec3_t *_P, const vec3_t *_Q, 
			   const mat33_t *F, const mat33_t *G, const int n);
void objpose(mat33_t *R, vec3_t *t, int *it, real_t *obj_err, real_t *img_err,
			 bool calc_img_err, const vec3_t *_P, const vec3_t *Qp, const options_t options, const int n);
void getRotationY_wrtT(double *al_ret, vec3_t *tnew, const vec3_t *v,
					   const vec3_t *p, const vec3_t *t, const real_t *DB,
					   const mat33_t *Rz, const int n, int *res_n);
void getRfor2ndPose_V_Exact(pose_t *sol, const vec3_t *v, const vec3_t *P,
					        const mat33_t R, const vec3_t t, const real_t DB, const int n, int *r_n);
void get2ndPose_Exact(pose_t *sol, const vec3_t *v, const vec3_t *P,
					  const mat33_t R, const vec3_t t, const real_t DB, const int n, int *sol_n);

// ------------------------------------------------------------------------------------------
void robust_pose(real_t *err, mat33_t *R , vec3_t *t, 
		 real_t *errb, mat33_t *Rb, vec3_t *tb,
		 const vec3_t *_model, const vec3_t *_iprts, const options_t _options, const int n);
void robust_pose2(real_t *err, mat33_t *R, vec3_t *t,
		 const vec3_t *_model, const vec3_t *_iprts,
		 const options_t _options, const int n);
// ------------------------------------------------------------------------------------------
//} // namespace rpp

#endif
