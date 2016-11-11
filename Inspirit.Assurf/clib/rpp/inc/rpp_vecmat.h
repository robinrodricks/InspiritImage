
#ifndef __RPP_VECMAT_H__
#define __RPP_VECMAT_H__

#include "rpp_const.h"
#include "rpp_types.h"

void _sin_cos(real_t a, real_t *s, real_t *c);
real_t _sin(real_t a);
real_t _cos(real_t a);
real_t _atan2(real_t a, real_t b);
real_t _abs(real_t a);
real_t _acos(real_t a);
real_t _sqrt(real_t a);
real_t _pow(real_t a, real_t b);

inline void mat33_assign(mat33_t *m,
				  const real_t m00, const real_t m01, const real_t m02,
				  const real_t m10, const real_t m11, const real_t m12,
				  const real_t m20, const real_t m21, const real_t m22);

inline void vec3_assign(vec3_t *v, const real_t x, const real_t y, const real_t z);
inline void vec3_clear(vec3_t *v);
inline void vec3_copy(vec3_t *a, const vec3_t *b);
inline void vec3_array_sum(vec3_t *v_sum2, const vec3_t *va, const int n);
inline void vec3_array_sum2(double *v_sum1, const vec3_t *va, const int n);
inline void vec3_array_pow2(vec3_t *va, const int n);
inline void vec3_div(vec3_t *va, const real_t n);
inline void vec3_div_vec(vec3_t *va, const vec3_t *vb);
inline void vec3_mult(vec3_t *va, const real_t n);
inline void vec3_mult_vec(vec3_t *va, const vec3_t *vb);
inline void vec3_add(vec3_t *va, const real_t f);
inline void vec3_add_vec(vec3_t *va, const vec3_t *vb);
inline void vec3_add_vec2(vec3_t *va, const vec3_t *vb, const vec3_t *vc);
inline void vec3_sub(vec3_t *va, const real_t f);
inline void vec3_sub_vec(vec3_t *va, const vec3_t *vb);
inline void vec3_sub_vec2(vec3_t *va, const vec3_t *vb, const vec3_t *vc);
inline real_t vec3_dot(const vec3_t *va, const vec3_t *vb);
inline void vec3_cross(vec3_t *va, const vec3_t *vb, const vec3_t *vc);
inline real_t vec3_norm(const vec3_t *v);
inline real_t vec3_sum(const vec3_t *v);
inline void vec3_array_add(vec3_t *va, const vec3_t *a, const int n);
inline void vec3_array_sub(vec3_t *va, const vec3_t *a, const int n);
inline void vec3_array_set(vec3_t *va, const vec3_t *a, const bool mask[3], const int n);
inline void vec3_array_mult(vec3_t *va, const double *c, const int n);
void vec3_array_mean(vec3_t *v_mean, const vec3_t *va, const int n);
inline void vec3_mul_vec3trans(mat33_t *m, const vec3_t *va, const vec3_t *vb);
inline real_t vec3trans_mul_vec3(const vec3_t *va, const vec3_t *vb);
inline void mat33_clear(mat33_t *m);
inline void mat33_copy(mat33_t *md, const mat33_t *ms);
inline void mat33_to_col_vec3(vec3_t *c0, vec3_t *c1, vec3_t *c2, const mat33_t *m);
inline void mat33_div(mat33_t *m, const real_t f);
inline void mat33_eye(mat33_t *m);
inline real_t mat33_sum(const mat33_t *m);
inline int mat33_all_zeros(const mat33_t *m);
inline void mat33_set_all_zeros(mat33_t *m);
void mat33_array_sum(mat33_t *s, const mat33_t *ma, const int n);
inline void mat33_sub_mat2(mat33_t *mr, const mat33_t *ma, const mat33_t *mb);
inline void mat33_sub(mat33_t *ma, const mat33_t *mb);
inline void mat33_add_mat2(mat33_t *mr, const mat33_t *ma, const mat33_t *mb);
inline void mat33_add(mat33_t *ma, const mat33_t *mb);
inline real_t mat33_det(const mat33_t *a);
inline void mat33_inv(mat33_t *mi, const mat33_t *ma);
inline void mat33_mult_mat2(mat33_t *m0, const mat33_t *m1, const mat33_t *m2);
inline void mat33_mult(mat33_t *mr, const real_t n);
inline void mat33_transpose(mat33_t *t, const mat33_t m);
inline void vec3_mult_mat(vec3_t *v0, const mat33_t *m1, const vec3_t *v2);
void vec3_array_mult_vec2(vec3_t *va, const mat33_t *m, const vec3_t *vb, const int n);
void mat33_svd2(mat33_t *u, mat33_t *s, mat33_t *v, const mat33_t *m);
void quat_mult(quat_t *q, const real_t s);
real_t quat_norm(const quat_t *q);
inline void mat33_from_quat(mat33_t *m, const quat_t *q);
void normRv(vec3_t *n, const vec3_t *v);
void normRv2(vec3_t *normR_v, const vec3_t *v, const int n);

int solveBiCubic(double *k1, double a, double b, double c, double d, double e);
int solve_polynomial(double *sol, const double *coefficients, const int n);
void scalar_array_pow(double *sa, const real_t f, const int n);
void scalar_array_negate(double *sa, const int n);
void scalar_array_assign(double *sa,
						 const	real_t f,
						 const int sz);
void scalar_array_add_vec(double *sa, const double *sb, const int n);
void scalar_array_clear(double *sa, const int n);
void scalar_array_atan2(double *sa, 
						const double *sb,
						const double *sc, const int n);

void scalar_array_div(double *sa, real_t f, const int n);
void scalar_array_div_vec(double *sa, const double *sb, const int n);
void scalar_array_mult(double *sa, real_t f, const int n);
void scalar_array_add(double *sa, real_t f, const int n);
void scalar_array_sub(double *sa, real_t f, const int n);
inline void mat33_pow2(mat33_t *m);

#endif
