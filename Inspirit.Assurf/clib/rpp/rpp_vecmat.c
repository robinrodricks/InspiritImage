#include "inc/rpp_vecmat.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

int svdcmp( double *a, int m,int n, double *w,double *v);
int quartic(double dd[5], double sol[4], double soli[4], int* Nsol);
int cubic(double A[4], double X[3], int* L);


void _sin_cos(real_t a, real_t *s, real_t *c)
{ 
	real_t angleTmp = a * 0.6366197723675814;

	int qi = (angleTmp + 0.5 - (real_t)(angleTmp < 0 ? 1 : 0)); // round(angle)
	angleTmp -= qi; // x=frac

	real_t xq = angleTmp * angleTmp;
	// less precision
	//                      var f0:Number = 1.0 + xq * (xq * 0.25360671639164339 - 1.2336979844380824);
	//                      var f1:Number = angleTmp * (1.5707963267948966 + xq * (  xq * 0.079679708649230657 - 0.64596348437163809));

	//more precision
	real_t f0 = 1.0 + xq * (xq * (0.25360671639164339 - 0.020427240364907607 * xq) - 1.2336979844380824);
	real_t f1 = angleTmp * (1.5707963267948966 + xq * (  xq * (0.079679708649230657 - 0.0046002309092153379 * xq) - 0.64596348437163809));

	real_t qi1 = qi & 1;
	real_t qi2 = qi & 2;
	real_t qi1_2 = 1.0 - qi2;
	qi2 = qi1 * qi1_2;
	qi1 = (1.0 - qi1) * qi1_2;

	*c = qi1 * f0 - qi2 * f1;
	*s = qi1 * f1 + qi2 * f0;
}

real_t _sin(real_t a)
{ 
	real_t angleTmp = a * 0.6366197723675814;

	int qi = (angleTmp + 0.5 - (real_t)(angleTmp < 0 ? 1 : 0)); // round(angle)
	angleTmp -= qi; // x=frac

	real_t xq = angleTmp * angleTmp;
	// less precision
	//                      var f0:Number = 1.0 + xq * (xq * 0.25360671639164339 - 1.2336979844380824);
	//                      var f1:Number = angleTmp * (1.5707963267948966 + xq * (  xq * 0.079679708649230657 - 0.64596348437163809));

	//more precision
	real_t f0 = 1.0 + xq * (xq * (0.25360671639164339 - 0.020427240364907607 * xq) - 1.2336979844380824);
	real_t f1 = angleTmp * (1.5707963267948966 + xq * (  xq * (0.079679708649230657 - 0.0046002309092153379 * xq) - 0.64596348437163809));

	real_t qi1 = qi & 1;
	real_t qi2 = qi & 2;
	real_t qi1_2 = 1.0 - qi2;
	qi2 = qi1 * qi1_2;
	qi1 = (1.0 - qi1) * qi1_2;

	//cos = qi1 * f0 - qi2 * f1;
	return qi1 * f1 + qi2 * f0;
}
real_t _cos(real_t a)
{ 
	real_t angleTmp = a * 0.6366197723675814;

	int qi = (angleTmp + 0.5 - (real_t)(angleTmp < 0 ? 1 : 0)); // round(angle)
	angleTmp -= qi; // x=frac

	real_t xq = angleTmp * angleTmp;
	// less precision
	//                      var f0:Number = 1.0 + xq * (xq * 0.25360671639164339 - 1.2336979844380824);
	//                      var f1:Number = angleTmp * (1.5707963267948966 + xq * (  xq * 0.079679708649230657 - 0.64596348437163809));

	//more precision
	real_t f0 = 1.0 + xq * (xq * (0.25360671639164339 - 0.020427240364907607 * xq) - 1.2336979844380824);
	real_t f1 = angleTmp * (1.5707963267948966 + xq * (  xq * (0.079679708649230657 - 0.0046002309092153379 * xq) - 0.64596348437163809));

	real_t qi1 = qi & 1;
	real_t qi2 = qi & 2;
	real_t qi1_2 = 1.0 - qi2;
	qi2 = qi1 * qi1_2;
	qi1 = (1.0 - qi1) * qi1_2;

	return qi1 * f0 - qi2 * f1;
	//return qi1 * f1 + qi2 * f0;
}

real_t _abs(real_t a) { return((a>0?a:-a)); }
//real_t _atan2(real_t a, real_t b) { return(::atan2(a,b)); }
real_t _atan2(real_t y, real_t x)
{
	double angle, r ;
	double const c3 = 0.1821 ;
	double const c1 = 0.9675 ;
	double abs_y = _abs(y) + 2.220446049250313e-16;

	if (x >= 0) {
		r = (x - abs_y) / (x + abs_y) ;
		angle = CONST_PI / 4 ;
	} else {
		r = (x + abs_y) / (abs_y - x) ;
		angle = 3 * CONST_PI / 4 ;
	}
	angle += (c3*r*r - c1) * r ;
	return (y < 0) ? - angle : angle ;
}
real_t _acos(real_t a) { return(acos(a)); }
//inline real_t _sqrt(real_t a) { return(::sqrt(a)); }
real_t _sqrt(real_t x) 
{ 
	if(x < 1e-8) return 0.0;
	
  /* 64-bit version */
  union {
    double x ;
    long long  i ;
  } u ;

  double xhalf = (double) 0.5 * x ;

  /* convert floating point value in RAW integer */
  u.x = x ;

  /* gives initial guess y0 */
	u.i = 0x5fe6ec85e7de30daLL - (u.i >> 1) ;

  /* two Newton steps */
  u.x = u.x * ( (double) 1.5  - xhalf*u.x*u.x) ;
  u.x = u.x * ( (double) 1.5  - xhalf*u.x*u.x) ;
  
  return x * u.x;
}
real_t _pow(real_t a, real_t b) { return(pow(a,b)); }

// ---------------------------------------------------------------------------

inline void mat33_assign(mat33_t *m,
				  const real_t m00, const real_t m01, const real_t m02,
				  const real_t m10, const real_t m11, const real_t m12,
				  const real_t m20, const real_t m21, const real_t m22)
{
	real_t *mm = &*m->m;
	mm[0] = m00; 	mm[1] = m01;	mm[2] = m02;
	mm[3] = m10; 	mm[4] = m11;	mm[5] = m12;
	mm[6] = m20; 	mm[7] = m21;	mm[8] = m22;
}

inline void vec3_assign(vec3_t *v, const real_t x, const real_t y, const real_t z)
{
	v->v[0] = x;
	v->v[1] = y;
	v->v[2] = z;
}

inline void vec3_clear(vec3_t *v)
{
	v->v[0] = 0;
	v->v[1] = 0;
	v->v[2] = 0;
}

inline void vec3_copy(vec3_t *a, const vec3_t *b)
{
	a->v[0] = b->v[0];
	a->v[1] = b->v[1];
	a->v[2] = b->v[2];
}

inline void vec3_array_sum(vec3_t *v_sum2, const vec3_t *va, const int n)
{
	vec3_clear(v_sum2);
	int i;
	real_t *vs = &*v_sum2->v;
	for(i = 0; i < n; i++)
	{
		const real_t *v = &*va[i].v;
		vs[0] += v[0];
		vs[1] += v[1];
		vs[2] += v[2];
	}
}

inline void vec3_array_sum2(real_t *v_sum1, const vec3_t *va, const int n)
{
	int i;
	for(i=0; i<n; i++)
	{
		const real_t *v = &*va[i].v;
		v_sum1[i] = (v[0] + v[1] + v[2]);
	}
}

inline void vec3_array_pow2(vec3_t *va, const int n)
{
	int i;
	for(i=0; i<n; i++)
	{
		real_t *v = &*va[i].v;
		v[0] *= v[0];
		v[1] *= v[1];
		v[2] *= v[2];
	}
}

inline void vec3_div(vec3_t *va, const real_t n)
{
	const real_t invn = (real_t)1.0 / n;
	va->v[0] *= invn;
	va->v[1] *= invn;
	va->v[2] *= invn;
}

inline void vec3_div_vec(vec3_t *va, const vec3_t *vb)
{
	va->v[0] /= vb->v[0];
	va->v[1] /= vb->v[1];
	va->v[2] /= vb->v[2];
}

inline void vec3_mult(vec3_t *va, const real_t n)
{
	va->v[0] *= n;
	va->v[1] *= n;
	va->v[2] *= n;
}

inline void vec3_mult_vec(vec3_t *va, const vec3_t *vb)
{
	va->v[0] *= vb->v[0];
	va->v[1] *= vb->v[1];
	va->v[2] *= vb->v[2];
}

inline void vec3_add(vec3_t *va, const real_t f)
{
	va->v[0] += f;
	va->v[1] += f;
	va->v[2] += f;
}

inline void vec3_add_vec(vec3_t *va, const vec3_t *vb)
{
	va->v[0] += vb->v[0];
	va->v[1] += vb->v[1];
	va->v[2] += vb->v[2];
}

inline void vec3_add_vec2(vec3_t *va, const vec3_t *vb, const vec3_t *vc)
{
	va->v[0] = vb->v[0] + vc->v[0];
	va->v[1] = vb->v[1] + vc->v[1];
	va->v[2] = vb->v[2] + vc->v[2];
}

inline void vec3_sub(vec3_t *va, const real_t f)
{
	va->v[0] -= f;
	va->v[1] -= f;
	va->v[2] -= f;
}

inline void vec3_sub_vec(vec3_t *va, const vec3_t *vb)
{
	va->v[0] -= vb->v[0];
	va->v[1] -= vb->v[1];
	va->v[2] -= vb->v[2];
}

inline void vec3_sub_vec2(vec3_t *va, const vec3_t *vb, const vec3_t *vc)
{
	va->v[0] = vb->v[0] - vc->v[0];
	va->v[1] = vb->v[1] - vc->v[1];
	va->v[2] = vb->v[2] - vc->v[2];
}

inline real_t vec3_dot(const vec3_t *va, const vec3_t *vb)
{
	return(va->v[0]*vb->v[0]+va->v[1]*vb->v[1]+va->v[2]*vb->v[2]);
}

inline void vec3_cross(vec3_t *va, const vec3_t *vb, const vec3_t *vc)
{
	va->v[0] = (vb->v[1] * vc->v[2] - vc->v[1] * vb->v[2]);
	va->v[1] = (vb->v[2] * vc->v[0] - vc->v[2] * vb->v[0]);
	va->v[2] = (vb->v[0] * vc->v[1] - vc->v[0] * vb->v[1]);
}

inline real_t vec3_norm(const vec3_t *v)
{
	return(_sqrt(v->v[0]*v->v[0] + v->v[1]*v->v[1] + v->v[2]*v->v[2]));
}

inline real_t vec3_sum(const vec3_t *v)
{
	return(v->v[0] + v->v[1] + v->v[2]);
}

inline void vec3_array_add(vec3_t *va, const vec3_t *a, const int n)
{
	int i;
	for(i = 0; i < n; i++)
	{
		va[i].v[0] += a->v[0];
		va[i].v[1] += a->v[1];
		va[i].v[2] += a->v[2];
	}
}

inline void vec3_array_sub(vec3_t *va, const vec3_t *a, const int n) 
{
	int i;
	for(i = 0; i < n; i++)
	{
		va[i].v[0] -= a->v[0];
		va[i].v[1] -= a->v[1];
		va[i].v[2] -= a->v[2];
	}
}

inline void vec3_array_set(vec3_t *va, const vec3_t *a, const int mask[3], const int n)
{
	int i;
	for(i = 0; i < n; i++)
	{
		if(mask[0] == 1) va[i].v[0] = a->v[0];
		if(mask[1] == 1) va[i].v[1] = a->v[1];
		if(mask[2] == 1) va[i].v[2] = a->v[2];
	}
}

inline void vec3_array_mult(vec3_t *va, const double *c, const int n)
{
	int i;
	for(i = 0; i < n; i++)
	{
		va[i].v[0] *= c[i];
		va[i].v[1] *= c[i];
		va[i].v[2] *= c[i];
	}
}

void vec3_array_mean(vec3_t *v_mean, const vec3_t *va, const int n)
{
	vec3_array_sum(v_mean, va, n);
	real_t l = (real_t) n;
	vec3_div(v_mean, l);
}


inline void vec3_mul_vec3trans(mat33_t *m, const vec3_t *va, const vec3_t *vb)
{
	real_t *mm = &*m->m;
	const real_t *vav = &*va->v;
	const real_t *vbv = &*vb->v;
	mm[0] = vav[0] * vbv[0];
	mm[1] = vav[0] * vbv[1];
	mm[2] = vav[0] * vbv[2];
	mm[3] = vav[1] * vbv[0];
	mm[4] = vav[1] * vbv[1];
	mm[5] = vav[1] * vbv[2];
	mm[6] = vav[2] * vbv[0];
	mm[7] = vav[2] * vbv[1];
	mm[8] = vav[2] * vbv[2];
}

inline real_t vec3trans_mul_vec3(const vec3_t *va, const vec3_t *vb)
{
	return(va->v[0] * vb->v[0] + va->v[1] * vb->v[1] + va->v[2] * vb->v[2]);
}


inline void mat33_clear(mat33_t *m)
{
	m->m[0] = 0;
	m->m[1] = 0;
	m->m[2] = 0;
	m->m[3] = 0;
	m->m[4] = 0;
	m->m[5] = 0;
	m->m[6] = 0;
	m->m[7] = 0;
	m->m[8] = 0;
}

inline void mat33_copy(mat33_t *md, const mat33_t *ms)
{
	/*md->m[0] = ms->m[0];
	md->m[1] = ms->m[1];
	md->m[2] = ms->m[2];
	md->m[3] = ms->m[3];
	md->m[4] = ms->m[4];
	md->m[5] = ms->m[5];
	md->m[6] = ms->m[6];
	md->m[7] = ms->m[7];
	md->m[8] = ms->m[8];*/
	memcpy(&*md->m, &*ms->m, 9 * sizeof(real_t));
}


inline void mat33_to_col_vec3(vec3_t *c0, vec3_t *c1, vec3_t *c2, const mat33_t *m)
{
	const real_t *mm = &*m->m;
	c0->v[0] = mm[0];
	c1->v[0] = mm[1];
	c2->v[0] = mm[2];
	
	c0->v[1] = mm[3];
	c1->v[1] = mm[4];
	c2->v[1] = mm[5];
	
	c0->v[2] = mm[6];
	c1->v[2] = mm[7];
	c2->v[2] = mm[8];
}


inline void mat33_div(mat33_t *m, const real_t f)
{
	const real_t invf = (real_t)1.0 / f;
	real_t *mm = &*m->m;
	mm[0] *= invf;
	mm[1] *= invf;
	mm[2] *= invf;
	mm[3] *= invf;
	mm[4] *= invf;
	mm[5] *= invf;
	mm[6] *= invf;
	mm[7] *= invf;
	mm[8] *= invf;
}

inline void mat33_eye(mat33_t *m)
{
	real_t *mm = &*m->m;
	mm[0] = 1;
	mm[1] = 0;
	mm[2] = 0;
	mm[3] = 0;
	mm[4] = 1;
	mm[5] = 0;
	mm[6] = 0;
	mm[7] = 0;
	mm[8] = 1;
}

inline real_t mat33_sum(const mat33_t *m)
{
	const real_t *mm = &*m->m;
	real_t sum = 0.0;
	sum += mm[0];
	sum += mm[1];
	sum += mm[2];
	sum += mm[3];
	sum += mm[4];
	sum += mm[5];
	sum += mm[6];
	sum += mm[7];
	sum += mm[8];
	
	return(sum);
}


inline int mat33_all_zeros(const mat33_t *m)
{
	const real_t *mm = &*m->m;
	if(mm[0] != 0) return 0;
	if(mm[1] != 0) return 0;
	if(mm[2] != 0) return 0;
	if(mm[3] != 0) return 0;
	if(mm[4] != 0) return 0;
	if(mm[5] != 0) return 0;
	if(mm[6] != 0) return 0;
	if(mm[7] != 0) return 0;
	if(mm[8] != 0) return 0;
	return 1;
}

inline void mat33_set_all_zeros(mat33_t *m)
{
	real_t *mm = &*m->m;
	mm[0] = 0;
	mm[1] = 0;
	mm[2] = 0;
	mm[3] = 0;
	mm[4] = 0;
	mm[5] = 0;
	mm[6] = 0;
	mm[7] = 0;
	mm[8] = 0;
}


void mat33_array_sum(mat33_t *s, const mat33_t *ma, const int n)
{
	mat33_clear(s);
	int i;
	real_t *sm = &*s->m;
	for(i = 0; i < n; i++)
	{
		const real_t *mm = &*ma[i].m;
		sm[0] += mm[0];
		sm[3] += mm[3];
		sm[6] += mm[6];
		
		sm[1] += mm[1];
		sm[4] += mm[4];
		sm[7] += mm[7];
		
		sm[2] += mm[2];
		sm[5] += mm[5];
		sm[8] += mm[8];
	}
}


inline void mat33_sub_mat2(mat33_t *mr, const mat33_t *ma, const mat33_t *mb)
{
	real_t *mrm = &*mr->m;
	const real_t *mam = &*ma->m;
	const real_t *mbm = &*mb->m;
	mrm[0] = mam[0] - mbm[0];
	mrm[1] = mam[1] - mbm[1];
	mrm[2] = mam[2] - mbm[2];
	mrm[3] = mam[3] - mbm[3];
	mrm[4] = mam[4] - mbm[4];
	mrm[5] = mam[5] - mbm[5];
	mrm[6] = mam[6] - mbm[6];
	mrm[7] = mam[7] - mbm[7];
	mrm[8] = mam[8] - mbm[8];
}

inline void mat33_sub(mat33_t *ma, const mat33_t *mb)
{	
	real_t *mam = &*ma->m;
	const real_t *mbm = &*mb->m;
	mam[0] -= mbm[0];
	mam[1] -= mbm[1];
	mam[2] -= mbm[2];
	mam[3] -= mbm[3];
	mam[4] -= mbm[4];
	mam[5] -= mbm[5];
	mam[6] -= mbm[6];
	mam[7] -= mbm[7];
	mam[8] -= mbm[8];
}

inline void mat33_add_mat2(mat33_t *mr, const mat33_t *ma, const mat33_t *mb)
{
	real_t *mrm = &*mr->m;
	const real_t *mam = &*ma->m;
	const real_t *mbm = &*mb->m;
	mrm[0] = mam[0] + mbm[0];
	mrm[1] = mam[1] + mbm[1];
	mrm[2] = mam[2] + mbm[2];
	mrm[3] = mam[3] + mbm[3];
	mrm[4] = mam[4] + mbm[4];
	mrm[5] = mam[5] + mbm[5];
	mrm[6] = mam[6] + mbm[6];
	mrm[7] = mam[7] + mbm[7];
	mrm[8] = mam[8] + mbm[8];
}

inline void mat33_add(mat33_t *ma, const mat33_t *mb)
{
	real_t *mam = &*ma->m;
	const real_t *mbm = &*mb->m;
	mam[0] += mbm[0];
	mam[1] += mbm[1];
	mam[2] += mbm[2];
	mam[3] += mbm[3];
	mam[4] += mbm[4];
	mam[5] += mbm[5];
	mam[6] += mbm[6];
	mam[7] += mbm[7];
	mam[8] += mbm[8];
}


inline real_t mat33_det(const mat33_t *a)
{
	real_t determinant = a->m[0]*a->m[4]*a->m[8] + a->m[1]*a->m[5]*a->m[6] +
							a->m[2]*a->m[3]*a->m[7] - a->m[6]*a->m[4]*a->m[2] - 
							a->m[7]*a->m[5]*a->m[0] - a->m[8]*a->m[3]*a->m[1];
	return(determinant);
}

inline void mat33_inv(mat33_t *mi, const mat33_t *ma)
{
	const real_t determinant = (real_t)1.0 / mat33_det(ma);
	
	const real_t ma_0 = ma->m[0];
	const real_t ma_1 = ma->m[1];
	const real_t ma_2 = ma->m[2];
	const real_t ma_3 = ma->m[3];
	const real_t ma_4 = ma->m[4];
	const real_t ma_5 = ma->m[5];
	const real_t ma_6 = ma->m[6];
	const real_t ma_7 = ma->m[7];
	const real_t ma_8 = ma->m[8];
	
	real_t *mim = &*mi->m;
	
	mim[0] = (ma_4*ma_8 - ma_5*ma_7)*determinant;
	mim[1] = (ma_2*ma_7 - ma_1*ma_8)*determinant;
	mim[2] = (ma_1*ma_5 - ma_2*ma_4)*determinant;

	mim[3] = (ma_5*ma_6 - ma_3*ma_8)*determinant;
	mim[4] = (ma_0*ma_8 - ma_2*ma_6)*determinant;
	mim[5] = (ma_2*ma_3 - ma_0*ma_5)*determinant;

	mim[6] = (ma_3*ma_7 - ma_4*ma_6)*determinant;
	mim[7] = (ma_1*ma_6 - ma_0*ma_7)*determinant;
	mim[8] = (ma_0*ma_4 - ma_1*ma_3)*determinant;
}

inline void mat33_mult_mat2(mat33_t *m0, const mat33_t *m1, const mat33_t *m2)
{
	const real_t m1_0 = m1->m[0];
	const real_t m1_1 = m1->m[1];
	const real_t m1_2 = m1->m[2];
	const real_t m1_3 = m1->m[3];
	const real_t m1_4 = m1->m[4];
	const real_t m1_5 = m1->m[5];
	const real_t m1_6 = m1->m[6];
	const real_t m1_7 = m1->m[7];
	const real_t m1_8 = m1->m[8];
	
	const real_t m2_0 = m2->m[0];
	const real_t m2_1 = m2->m[1];
	const real_t m2_2 = m2->m[2];
	const real_t m2_3 = m2->m[3];
	const real_t m2_4 = m2->m[4];
	const real_t m2_5 = m2->m[5];
	const real_t m2_6 = m2->m[6];
	const real_t m2_7 = m2->m[7];
	const real_t m2_8 = m2->m[8];
	
	real_t *m0m = &*m0->m;
	
	m0m[0] = m1_0*m2_0 + m1_1*m2_3 + m1_2*m2_6;
	m0m[1] = m1_0*m2_1 + m1_1*m2_4 + m1_2*m2_7;
	m0m[2] = m1_0*m2_2 + m1_1*m2_5 + m1_2*m2_8;
	m0m[3] = m1_3*m2_0 + m1_4*m2_3 + m1_5*m2_6;
	m0m[4] = m1_3*m2_1 + m1_4*m2_4 + m1_5*m2_7;
	m0m[5] = m1_3*m2_2 + m1_4*m2_5 + m1_5*m2_8;
	m0m[6] = m1_6*m2_0 + m1_7*m2_3 + m1_8*m2_6;
	m0m[7] = m1_6*m2_1 + m1_7*m2_4 + m1_8*m2_7;
	m0m[8] = m1_6*m2_2 + m1_7*m2_5 + m1_8*m2_8;
}

inline void mat33_mult(mat33_t *mr, const real_t n)
{
	real_t *mrm = &*mr->m;
	mrm[0] *= n;
	mrm[1] *= n;
	mrm[2] *= n;
	mrm[3] *= n;
	mrm[4] *= n;
	mrm[5] *= n;
	mrm[6] *= n;
	mrm[7] *= n;
	mrm[8] *= n;
}

inline void mat33_transpose(mat33_t *t, const mat33_t m)
{	
	const real_t *mm = &*m.m;
	real_t *tm = &*t->m;
	tm[0] = mm[0];
	tm[1] = mm[3];
	tm[2] = mm[6];
	
	tm[3] = mm[1];
	tm[4] = mm[4];
	tm[5] = mm[7];
	
	tm[6] = mm[2];
	tm[7] = mm[5];
	tm[8] = mm[8];
}

inline void vec3_mult_mat(vec3_t *v0, const mat33_t *m1, const vec3_t *v2)
{
	v0->v[0] = m1->m[0]*v2->v[0] + m1->m[1]*v2->v[1] + m1->m[2]*v2->v[2];
	v0->v[1] = m1->m[3]*v2->v[0] + m1->m[4]*v2->v[1] + m1->m[5]*v2->v[2];
	v0->v[2] = m1->m[6]*v2->v[0] + m1->m[7]*v2->v[1] + m1->m[8]*v2->v[2];
}

void vec3_array_mult_vec2(vec3_t *va, const mat33_t *m, const vec3_t *vb, const int n)
{
	int i;
	for(i=0; i<n; i++)
	{
		vec3_mult_mat(&va[i], m, &vb[i]);
	}
}

void mat33_svd2(mat33_t *u, mat33_t *s, mat33_t *v, const mat33_t *m)
{
	mat33_clear(u);
	mat33_clear(v);

	real_t m_ptr[9];
	
	double q_ptr[3] = {0.0, 0.0, 0.0};
	
	memcpy(&*m_ptr, &*m->m, 9 * sizeof(real_t));
	svdcmp(&*m_ptr, 3, 3, &*q_ptr, &*v->m);
	
	memcpy(&*u->m, &*m_ptr, 9 * sizeof(real_t));

	mat33_clear(s);
	real_t *sm = &*s->m;
	real_t *um = &*u->m;
	real_t *vm = &*v->m;
	sm[0] = (real_t)q_ptr[0];
	sm[4] = (real_t)q_ptr[1];
	sm[8] = (real_t)q_ptr[2];
	
	// WE need to sort the diagonal values of the result
	// this is neccessery for the absolute orientation algo
	// Biggest value first
	
	int sorted = 0;
	real_t t;
	while( sorted == 0 )
	{
	  
		sorted = 1;

		if( sm[0] < sm[4] )
		{
			sorted = 0;   
			// change i mit i-1 !!!

			// in the S
			t = sm[0];
			sm[0]  =  sm[4];
			sm[4] = t;
			
			// in U 
			t = um[0];  
			um[0] =  um[1];  
			um[1] = t;
			// and in V
			t = vm[0];  
			vm[0] = vm[1];  
			vm[1] = t;
			
			// in U 
			t = um[3];  
			um[3] =  um[4];  
			um[4] = t;
			// and in V
			t = vm[3];  
			vm[3] = vm[4];  
			vm[4] = t;
			
			// in U 
			t = um[6];  
			um[6] =  um[7];  
			um[7] = t;
			// and in V
			t = vm[6];  
			vm[6] = vm[7];  
			vm[7] = t;
		}
		//
		if( sm[4] < sm[8] )
		{
			sorted = 0;   
			// change i mit i-1 !!!

			// in the S
			t =  sm[4];
			sm[4] = sm[8];
			sm[8] = t;

			// in U 
			t = um[1];  
			um[1] = um[2];  
			um[2] = t;
			// and in V
			t = vm[1];  
			vm[1] = vm[2];  
			vm[2] = t;
			
			// in U 
			t = um[4];  
			um[4] = um[5];  
			um[5] = t;
			// and in V
			t = vm[4];  
			vm[4] = vm[5];  
			vm[5] = t;
			
			// in U 
			t = um[7];  
			um[7] = um[8];  
			um[8] = t;
			// and in V
			t = vm[7];  
			vm[7] = vm[8];  
			vm[8] = t;
		}
		
	}
}

void quat_mult(quat_t *q, const real_t s)
{
	vec3_mult(&q->v, s);
	q->s *= s;
}

real_t quat_norm(const quat_t *q)
{
	const real_t f_vn = vec3_norm(&q->v);
	return(_sqrt((f_vn*f_vn) + (q->s*q->s)));
}

inline void mat33_from_quat(mat33_t *m, const quat_t *q)
{
	const real_t a = q->s;
	const real_t b = q->v.v[0];
	const real_t c = q->v.v[1];
	const real_t d = q->v.v[2];
	real_t *mm = &*m->m;

	mm[0] = (a*a)+(b*b)-(c*c)-(d*d);
	mm[1] = (real_t)(2.0)*(b*c-a*d);
	mm[2] = (real_t)(2.0)*(b*d+a*c);

	mm[3] = (real_t)(2.0)*(b*c+a*d);
	mm[4] = (a*a)+(c*c)-(b*b)-(d*d);
	mm[5] = (real_t)(2.0)*(c*d-a*b);

	mm[6] = (real_t)(2.0)*(b*d-a*c);
	mm[7] = (real_t)(2.0)*(c*d+a*b);
	mm[8] = (a*a)+(d*d)-(b*b)-(c*c);
}

// ===========================================================================================
void normRv(vec3_t *n, const vec3_t *v)
{
	vec3_t _v1;
	vec3_copy(&_v1, v);
	vec3_mult_vec(&_v1, &_v1);
	real_t l = (real_t)1.0 / _sqrt(_v1.v[0] + _v1.v[1] + _v1.v[2]);
	vec3_copy(n, v);
	vec3_mult(n, l);
}

// ===========================================================================================
void normRv2(vec3_t *normR_v, const vec3_t *v, const int n)
{
	//normR_v.assign(v.begin(),v.end());
	memcpy(normR_v, v, n*sizeof(vec3_t));
	vec3_array_pow2(normR_v, n);
	double _l[n];
	vec3_array_sum2(&*_l, normR_v, n);
	
	int i;
	for(i = 0; i < n; i++)
	{
		_l[i] = (real_t)1.0 / _sqrt(_l[i]);
	}

	//normR_v.assign(v.begin(),v.end());
	memcpy(normR_v, v, n*sizeof(vec3_t));
	vec3_array_mult(normR_v, &*_l, n);
}

// ===========================================================================================
int solve_polynomial(double *r_sol, const double *coefficients, const int n)
{
	double dd[5] = {(double)coefficients[0],
		(double)coefficients[1],
		(double)coefficients[2],
		(double)coefficients[3],
		(double)coefficients[4] };

	double sol[4] = {0,0,0,0};
	double soli[4] = {0,0,0,0};
	int n_sol = 0;
	quartic(&*dd, &*sol, &*soli, &n_sol);

	if(n_sol <= 0) return(0);

	//int i;
	//for( i=0; i<n_sol; i++) r_sol[i] = (real_t)sol[i];
	memcpy(r_sol, &*sol, n_sol * sizeof(real_t));
	return(n_sol);
}
// ===========================================================================================

void scalar_array_pow(double *sa, const real_t f, const int n)
{
	int i;
	for( i=0; i<n; i++) sa[i] = _pow(sa[i], f);
}

void scalar_array_negate(double *sa, const int n)
{
	int i;
	for(i=0; i<n; i++) sa[i] = - sa[i];
}

void scalar_array_assign(double *sa,
						 const	real_t f,
						 const int sz)
{
	int i;
	for( i=0; i<sz; i++) sa[i] = f;
}


void scalar_array_add_vec(double *sa, const double *sb, const int n)
{
	int i;
	for(i=0; i<n; i++)
	{
		sa[i] = sa[i] + sb[i];
	}
}

void scalar_array_clear(double *sa, const int n)
{
	int i;
	for(i=0; i<n; i++) sa[i] = 0.;
}

void scalar_array_atan2(double *sa, 
						const double *sb,
						const double *sc, const int n)
{
	int i;
	for(i=0; i<n; i++)
	{
		sa[i] = _atan2(sb[i], sc[i]);
	}
}

void scalar_array_div(double *sa, real_t f, const int n)
{
	int i;
	const real_t invf = (real_t)1.0 / f;
	for(i=0; i<n; i++)
	{
		sa[i] *= invf;
	}
}

void scalar_array_div_vec(double *sa, const double *sb, const int n)
{
	int i;
	for(i=0; i<n; i++)
	{
		sa[i] /= sb[i];
	}
}

void scalar_array_mult(double *sa, real_t f, const int n)
{
	int i;
	for(i=0; i<n; i++)
	{
		sa[i] *= f;
	}
}

void scalar_array_add(double *sa, real_t f, const int n)
{
	int i;
	for(i=0; i<n; i++)
	{
		sa[i] += f;
	}
}

void scalar_array_sub(double *sa, real_t f, const int n)
{
	int i;
	for(i=0; i<n; i++)
	{
		sa[i] -= f;
	}
}

inline void mat33_pow2(mat33_t *m)
{
	real_t *mm = &*m->m;
	mm[0] *= mm[0];
	mm[1] *= mm[1];
	mm[2] *= mm[2];
	mm[3] *= mm[3];
	mm[4] *= mm[4];
	mm[5] *= mm[5];
	mm[6] *= mm[6];
	mm[7] *= mm[7];
	mm[8] *= mm[8];
}
