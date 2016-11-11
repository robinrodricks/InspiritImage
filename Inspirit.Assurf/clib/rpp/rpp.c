#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "inc/rpp.h"
#include "inc/rpp_const.h"
#include "inc/rpp_vecmat.h"

#define CALC_IMAGE_ERROR 0
#define COMPLICATED_ERROR_CMP 0

// ===========================================================================================
void Quaternion_byAngleAndVector(quat_t *Q, const real_t q_angle, const vec3_t *q_vector)
{
	vec3_t rotation_axis;
	normRv(&rotation_axis, q_vector);
	
	//real_t f = _sin(q_angle * 0.5);
	real_t f, cs;
	_sin_cos(q_angle*0.5, &f, &cs);
	
	vec3_copy(&Q->v, &rotation_axis);
	vec3_mult(&Q->v, f);
	//Q->s = _cos(q_angle*0.5);
	Q->s = cs;
	quat_mult(Q, 1.0 / quat_norm(Q));
}

// ===========================================================================================

void GetRotationbyVector(mat33_t *R, const vec3_t *v1, const vec3_t *v2)
{
	 vec3_t diff;
  vec3_sub_vec2(&diff, v1, v2);
    
  if( vec3_norm(&diff) < 1e-3 ){
    mat33_assign(R,1,0,0,0,1,0,0,0,1);
  }else
  {
		if( fabs((vec3_norm(&diff)-2)) < 1e-3 ){
		  mat33_assign(R,1,0,0 ,0,-1,0, 0,0,-1);
		}else
		  { 
			real_t winkel = _acos(vec3_dot(v1,v2));
			quat_t QU;
			vec3_t vc;
			vec3_cross(&vc,v2,v1);
			Quaternion_byAngleAndVector(&QU,winkel,&vc);
			mat33_from_quat(R,&QU);
		}
	}
}

// ===========================================================================================
void xform(vec3_t* Q, const vec3_t *P, const mat33_t *R, const vec3_t *t, const int n)
{
	//const unsigned int n = (unsigned int) P.size();
	int i;
	for(i=0; i<n; i++)
	{
		vec3_mult_mat(&Q[i], R, &P[i]);
		vec3_add_vec(&Q[i], t);
	}
}
// ===========================================================================================

void xformproj(vec3_t *Qp, const vec3_t *P, const mat33_t *R, const vec3_t *t, const int n)
{
	//const unsigned int n = (unsigned int) P.size();
	vec3_t Q[n];
	//Q.resize(n);
	int i;
	for(i=0; i<n; i++)
	{
		vec3_mult_mat(&Q[i], R, &P[i]);
		vec3_add_vec(&Q[i], t);
		Qp[i].v[0] = Q[i].v[0] / Q[i].v[2]; 
		Qp[i].v[1] = Q[i].v[1] / Q[i].v[2]; 
		Qp[i].v[2] = 1.0;
	}
}

// ===========================================================================================
void rpyMat(mat33_t *R, const vec3_t *rpy) // rpy: roll,pitch,yaw
{
	/*const real_t cosA = _cos(rpy->v[2]);
	const real_t sinA = _sin(rpy->v[2]);
	const real_t cosB = _cos(rpy->v[1]);
	const real_t sinB = _sin(rpy->v[1]);
	const real_t cosC = _cos(rpy->v[0]);
	const real_t sinC = _sin(rpy->v[0]);*/
	
	real_t cosA, sinA, cosB, sinB, cosC, sinC;
	_sin_cos(rpy->v[2], &sinA, &cosA);
	_sin_cos(rpy->v[1], &sinB, &cosB);
	_sin_cos(rpy->v[0], &sinC, &cosC);
	
	
	const real_t cosAsinB = cosA * sinB;
	const real_t sinAsinB = sinA * sinB;

	R->m[0] = cosA*cosB;
	R->m[1] = cosAsinB*sinC-sinA*cosC;
	R->m[2] = cosAsinB*cosC+sinA*sinC;

	R->m[3] = sinA*cosB;
	R->m[4] = sinAsinB*sinC+cosA*cosC;
	R->m[5] = sinAsinB*cosC-cosA*sinC;

	R->m[6] = -sinB;
	R->m[7] = cosB*sinC;
	R->m[8] = cosB*cosC;
}
// ===========================================================================================
void rpyAng(vec3_t *angs, const mat33_t *R)
{
	const real_t sinB = -(R->m[6]);
	const real_t cosB = _sqrt(R->m[0]*R->m[0] + R->m[3]*R->m[3]);

	if(_abs(cosB) > (real_t)(1E-15))
	{
		const real_t sinA = R->m[3] / cosB;
		const real_t cosA = R->m[0] / cosB;
		const real_t sinC = R->m[7] / cosB;
		const real_t cosC = R->m[8] / cosB;
		vec3_assign(angs,_atan2(sinC,cosC),_atan2(sinB,cosB),_atan2(sinA,cosA));
	}
	else
	{
		const real_t sinC = (R->m[1] - R->m[5]) * 0.5;
		const real_t cosC = (R->m[4] - R->m[2]) * 0.5;
		vec3_assign(angs,_atan2(sinC,cosC),CONST_PI_OVER_2, 0.0);
	}
}

// ===========================================================================================

void rpyAng_X(vec3_t *ang_zyx, const mat33_t *R)
{
	rpyAng(ang_zyx,R);

	if(_abs(ang_zyx->v[0]) > CONST_PI_OVER_2)
	{
		while(_abs(ang_zyx->v[0]) > CONST_PI_OVER_2)
		{
			if(ang_zyx->v[0] > 0)
			{
				vec3_assign(ang_zyx, ang_zyx->v[0]+CONST_PI,
					                 3*CONST_PI-ang_zyx->v[1],
									 ang_zyx->v[2]+CONST_PI);
				vec3_sub(ang_zyx,CONST_2_PI);
			}
			else
			{
				vec3_assign(ang_zyx, ang_zyx->v[0]+CONST_PI,
									 3*CONST_PI-ang_zyx->v[1],
									 ang_zyx->v[2]+CONST_PI);
			}
		}
	}
}

// ===========================================================================================
void decomposeR(mat33_t *Rz, const mat33_t *R)
{
	real_t cl = _atan2(R->m[7],R->m[6]);
	vec3_t rpy;
	vec3_assign(&rpy,0,0,cl);
	rpyMat(Rz,&rpy);
}


void optimal_t(vec3_t *t ,const mat33_t *R, const mat33_t *G,
	       const mat33_t *F,const vec3_t *P, int n){
  
  vec3_t _sum;
  vec3_clear(&_sum);
  
  int i;
  for(i=0; i<n; i++)
  {
    vec3_t _v1, _v2;
    vec3_mult_mat(&_v1, R, &P[i]);
    vec3_mult_mat(&_v2, &F[i], &_v1);
    vec3_add_vec(&_sum, &_v2);
  }
  
  vec3_mult_mat(t, G, &_sum);
}


// ===========================================================================================
void abskernel(mat33_t *R, vec3_t *t, vec3_t *Qout, real_t *err2, 
			   const vec3_t *_P, const vec3_t *_Q, 
			   const mat33_t *F, const mat33_t *G, const int n)

{
	int i, j;

	vec3_t P[n];//(_P.begin(),_P.end());
	vec3_t Q[n];//(_Q.begin(),_Q.end());
	//const unsigned int n = (unsigned int) P.size();
	memcpy(&*P, _P, n*sizeof(vec3_t));
	memcpy(&*Q, _Q, n*sizeof(vec3_t));

	for(i=0; i<n; i++)
	{
		vec3_t _q;
		vec3_mult_mat(&_q, &F[i], &_Q[i]);
		vec3_copy(&Q[i], &_q);
	}

	vec3_t pbar;
	vec3_array_sum(&pbar, &*P, n);
	vec3_div(&pbar, (real_t)(n));
	vec3_array_sub(&*P, &pbar, n);

	vec3_t qbar;
	vec3_array_sum(&qbar, &*Q, n);
	vec3_div(&qbar, (real_t)(n));
	vec3_array_sub(&*Q, &qbar, n);

	mat33_t M;
	mat33_clear(&M);
	for(j=0; j<n; j++)
	{
		mat33_t _m;
		vec3_mul_vec3trans(&_m, &P[j], &Q[j]);
		mat33_add(&M, &_m);
	}

	mat33_t _U;
	mat33_t _S;
	mat33_t _V;
	mat33_clear(&_U);
	mat33_clear(&_S);
	mat33_clear(&_V);
	mat33_svd2(&_U, &_S, &_V, &M);

	mat33_t _Ut;
	mat33_transpose(&_Ut, _U);
	mat33_mult_mat2(R, &_V, &_Ut);

	// R = V*(U.');

	// WE need to check the determinant of the R
	real_t det = mat33_det(R);
	
	
	if ( det < 0 )
	{
	  mat33_t Vc; 
	  Vc.m[0] = _V.m[0]; Vc.m[1] = _V.m[1]; Vc.m[2] = -_V.m[2];
	  Vc.m[3] = _V.m[3]; Vc.m[4] = _V.m[4]; Vc.m[5] = -_V.m[5];
	  Vc.m[6] = _V.m[6]; Vc.m[7] = _V.m[7]; Vc.m[8] = -_V.m[8];
	  mat33_mult_mat2(R, &Vc, &_Ut);  // to have det(R) == 1
	  optimal_t(t, R, G, F, &*P, n);
	  if( t->v[2] < 0 ){
	    Vc.m[0] =- _V.m[0]; Vc.m[1] =- _V.m[1]; Vc.m[2] = -_V.m[2];
	    Vc.m[3] =- _V.m[3]; Vc.m[4] =- _V.m[4]; Vc.m[5] = -_V.m[5];
	    Vc.m[6] =- _V.m[6]; Vc.m[7] =- _V.m[7]; Vc.m[8] = -_V.m[8];
	    mat33_mult_mat2(R, &Vc, &_Ut);  // to have det(R) == 1 & t_3 > 0
	    optimal_t(t, R, G, F, &*P, n);
	  }
	}else{
	  
	  optimal_t(t,R,G,F,&*P,n);
	  if( t->v[2] < 0 )
	  {
	    
	    mat33_t Vc; 
	    Vc.m[0] = -_V.m[0]; Vc.m[1] = -_V.m[1]; Vc.m[2] = _V.m[2];
        Vc.m[3] = -_V.m[3]; Vc.m[4] = -_V.m[4]; Vc.m[5] = _V.m[5];
        Vc.m[6] = -_V.m[6]; Vc.m[7] = -_V.m[7]; Vc.m[8] = _V.m[8];
	    mat33_mult_mat2(R, &Vc, &_Ut);
	    optimal_t(t,R,G,F,&*P,n);
	    
	  }
	}
	
	// CHECK if everything is working ok 
	// check the det of R 
	//det = mat33_det(R);
	
	xform(Qout,&*P,R,t, n);
	*err2 = 0.0;
	mat33_t _m1;
	vec3_t _v1;
	for(i=0; i<n; i++)
	{
		mat33_eye(&_m1);
		mat33_sub(&_m1, &F[i]);
		vec3_mult_mat(&_v1, &_m1, &Qout[i]);
		*err2 += vec3_dot(&_v1, &_v1);
	}
}
// ===========================================================================================

void objpose(mat33_t *R, vec3_t *t, int *it, real_t *obj_err, real_t *img_err,
			 bool calc_img_err, const vec3_t *_P, const vec3_t *Qp, const options_t options, const int n)
{
	int i, j;
	//vec3_array P(_P.begin(),_P.end());
	vec3_t P[n];
	memcpy(&*P, _P, n*sizeof(vec3_t));

	//const int n = (unsigned int) P.size();
	vec3_t pbar;
	vec3_array_sum(&pbar, &*P, n);
	vec3_div(&pbar, (real_t)(n));
	vec3_array_sub(&*P, &pbar, n);
	
	//vec3_array Q(Qp.begin(),Qp.end());
	vec3_t Q[n];
	memcpy(&*Q, Qp, n*sizeof(vec3_t));
	
	vec3_t ones;
	ones.v[0] = 1;
	ones.v[1] = 1;
	ones.v[2] = 1;
	const bool mask_z[3] = {0,0,1};
	vec3_array_set(&*Q, &ones, mask_z, n);
	
	//mat33_array F;
	//F.resize(n);
	mat33_t F[n];
	
	vec3_t V;
	for(i=0; i<n; i++)
	{
		V.v[0] = Q[i].v[0] / Q[i].v[2];
		V.v[1] = Q[i].v[1] / Q[i].v[2];
		V.v[2] = 1.0;
		mat33_t _m;
		vec3_mul_vec3trans(&_m, &V, &V);
		mat33_div(&_m, vec3trans_mul_vec3(&V,&V));
		F[i] = _m;
	}

	mat33_t tFactor;
	
	mat33_t _m1,_m2,_m3;
	mat33_eye(&_m1);
	mat33_array_sum(&_m2, &*F, n);
	mat33_div(&_m2, (real_t)(n));
	mat33_sub_mat2(&_m3, &_m1, &_m2);
	mat33_inv(&tFactor, &_m3);
	mat33_div(&tFactor, (real_t)(n));

	*it = 0;
	int initR_approximate = mat33_all_zeros(&options.initR);
	mat33_t Ri;
	vec3_t ti;
	
	//vec3_array Qi;
	//Qi.resize(n);
	vec3_t Qi[n];
	
	real_t old_err = 0.0, new_err = 0.0;

	// ----------------------------------------------------------------------------------------
	if(initR_approximate == 0)
	{
		mat33_copy(&Ri, &options.initR);
		vec3_t _sum;
		vec3_t _v1, _v2;
		mat33_t _m1,_m2;
		vec3_clear(&_sum);
		for(j=0; j<n; j++)
		{
			mat33_eye(&_m1);
			mat33_sub_mat2(&_m2, &F[j], &_m1);
			vec3_mult_mat(&_v1, &Ri, &P[j]);
			vec3_mult_mat(&_v2, &_m2, &_v1);
			vec3_add_vec(&_sum, &_v2);
		}
		vec3_mult_mat(&ti,&tFactor,&_sum);
		xform(&*Qi, &*P, &Ri, &ti, n);
		old_err = 0;
		vec3_t _v;
		for(j=0; j<n; j++)
		{
			mat33_eye(&_m1);
			mat33_sub_mat2(&_m2, &F[j], &_m1);
			vec3_mult_mat(&_v, &_m2, &Qi[j]);
			old_err += vec3_dot(&_v, &_v);
		}
	// ----------------------------------------------------------------------------------------
	}
	else
	{
		abskernel(&Ri, &ti, &*Qi, &old_err, &*P, &*Q, &*F, &tFactor, n);
		*it = 1;
	}
	// ----------------------------------------------------------------------------------------

	abskernel(&Ri, &ti, &*Qi, &new_err, &*P, &*Qi, &*F, &tFactor, n);
	*it = *it + 1;

	while((_abs((old_err-new_err)/old_err) > options.tol) && (new_err > options.epsilon) &&
		  (options.max_iter == 0 || *it < options.max_iter))
	{
		old_err = new_err;
		abskernel(&Ri, &ti, &*Qi, &new_err, &*P, &*Qi, &*F, &tFactor, n);
		*it = *it + 1;
	}


	mat33_copy(R, &Ri);
	vec3_copy(t, &ti);
	*obj_err = _sqrt(new_err/(real_t)(n));

	if(calc_img_err == 1)
	{
		//vec3_array Qproj;
		//Qproj.resize(n);
		vec3_t Qproj[n];
		
		xformproj(&*Qproj, &*P, &Ri, &ti, n);
		*img_err = 0;

		vec3_t _v;
		for(j=0; j<n; j++)
		{
			vec3_sub_vec2(&_v, &Qproj[j], &Qp[j]);
			*img_err += vec3_dot(&_v, &_v);
		}
		*img_err = _sqrt(*img_err/(real_t)(n));
	}

	if(t->v[2] < 0)
	{
		mat33_mult(R, -1.0);
		vec3_mult(t, -1.0);
	}

	vec3_t _ts;
	vec3_mult_mat(&_ts, &Ri, &pbar);
	vec3_sub_vec(t, &_ts);
}

// =====================================================================================

void getRotationY_wrtT(double *al_ret, vec3_t *tnew, const vec3_t *v,
					   const vec3_t *p, const vec3_t *t, const real_t *DB,
					   const mat33_t *Rz, const int n, int *res_n)
{
	int i,j,k;
	
	mat33_t V[n];
	
	for(i=0; i<n; i++)
	{
		vec3_mul_vec3trans(&V[i], &v[i], &v[i]);
		mat33_div(&V[i], vec3trans_mul_vec3(&v[i], &v[i]));
	}

	mat33_t G, _g1, _g2, _g3;
	mat33_array_sum(&_g1, &*V, n);
	mat33_eye(&_g2);
	mat33_div(&_g1, (real_t)(n));
	mat33_sub_mat2(&_g3, &_g2, &_g1);
	mat33_inv(&G, &_g3);
	mat33_div(&G, (real_t)(n));
	mat33_t _opt_t;
	mat33_clear(&_opt_t);

	for(i=0; i<n; i++)
	{
		const real_t v11 = V[i].m[0]; 
		const real_t v21 = V[i].m[3];
		const real_t v31 = V[i].m[6];
		const real_t v12 = V[i].m[1]; 
		const real_t v22 = V[i].m[4];
		const real_t v32 = V[i].m[7];
		const real_t v13 = V[i].m[2]; 
		const real_t v23 = V[i].m[5];
		const real_t v33 = V[i].m[8];
		const real_t px = p[i].v[0];
		const real_t py = p[i].v[1];
		const real_t pz = p[i].v[2];
		const real_t r1 = Rz->m[0];
		const real_t r2 = Rz->m[1];
		const real_t r3 = Rz->m[2];
		const real_t r4 = Rz->m[3];
		const real_t r5 = Rz->m[4];
		const real_t r6 = Rz->m[5];
		const real_t r7 = Rz->m[6];
		const real_t r8 = Rz->m[7];
		const real_t r9 = Rz->m[8];

		mat33_t _o;
		_o.m[0] = (((v11-(real_t)(1))*r2+v12*r5+v13*r8)*py+(-(v11-(real_t)(1))*r1-v12*r4-v13*r7)*px+(-(v11-(real_t)(1))*r3-v12*r6-v13*r9)*pz);
		_o.m[1] = (((real_t)(2)*(v11-(real_t)(1))*r1+(real_t)(2)*v12*r4+(real_t)(2)*v13*r7)*pz+(-(real_t)(2)*(v11-(real_t)(1))*r3-(real_t)(2)*v12*r6-(real_t)(2)*v13*r9)*px);
		_o.m[2] = ((v11-(real_t)(1))*r1+v12*r4+v13*r7)*px+((v11-(real_t)(1))*r3+v12*r6+v13*r9)*pz+((v11-(real_t)(1))*r2+v12*r5+v13*r8)*py;

		_o.m[3] = ((v21*r2+(v22-(real_t)(1))*r5+v23*r8)*py+(-v21*r1-(v22-(real_t)(1))*r4-v23*r7)*px+(-v21*r3-(v22-(real_t)(1))*r6-v23*r9)*pz);
		_o.m[4] = (((real_t)(2)*v21*r1+(real_t)(2)*(v22-(real_t)(1))*r4+(real_t)(2)*v23*r7)*pz+(-(real_t)(2)*v21*r3-(real_t)(2)*(v22-(real_t)(1))*r6-(real_t)(2)*v23*r9)*px);
		_o.m[5] = (v21*r1+(v22-(real_t)(1))*r4+v23*r7)*px+(v21*r3+(v22-(real_t)(1))*r6+v23*r9)*pz+(v21*r2+(v22-(real_t)(1))*r5+v23*r8)*py;

		_o.m[6] = ((v31*r2+v32*r5+(v33-(real_t)(1))*r8)*py+(-v31*r1-v32*r4-(v33-(real_t)(1))*r7)*px+(-v31*r3-v32*r6-(v33-(real_t)(1))*r9)*pz);
		_o.m[7] = (((real_t)(2)*v31*r1+(real_t)(2)*v32*r4+(real_t)(2)*(v33-(real_t)(1))*r7)*pz+(-(real_t)(2)*v31*r3-(real_t)(2)*v32*r6-(real_t)(2)*(v33-(real_t)(1))*r9)*px);
		_o.m[8] = (v31*r1+v32*r4+(v33-(real_t)(1))*r7)*px+(v31*r3+v32*r6+(v33-(real_t)(1))*r9)*pz+(v31*r2+v32*r5+(v33-(real_t)(1))*r8)*py;

		mat33_add(&_opt_t, &_o);
	}

	mat33_t opt_t;
	mat33_mult_mat2(&opt_t, &G, &_opt_t);
	real_t E_2[5] = {0,0,0,0,0};
	for(i=0; i<n; i++)
	{
		const real_t px = p[i].v[0];
		const real_t py = p[i].v[1];
		const real_t pz = p[i].v[2];

		mat33_t Rpi;
		mat33_assign(&Rpi, -px, (real_t)(2)*pz,px,py,(real_t)(0),py,-pz,-(real_t)(2)*px,pz);

		mat33_t E,_e1,_e2;
		mat33_eye(&_e1);
		mat33_sub(&_e1, &V[i]);
		mat33_mult_mat2(&_e2, Rz, &Rpi);
		mat33_add(&_e2, &opt_t);
		mat33_mult_mat2(&E,&_e1,&_e2);
		vec3_t e2,e1,e0;
		mat33_to_col_vec3(&e2,&e1,&e0,&E);
		vec3_t _E2_0,_E2_1,_E2_2,_E2_3,_E2_4;
		vec3_copy(&_E2_0,&e2);
		vec3_mult_vec(&_E2_0,&e2);
		vec3_copy(&_E2_1,&e1);
		vec3_mult_vec(&_E2_1,&e2);
		vec3_mult(&_E2_1,2.0);
		vec3_copy(&_E2_2,&e0);
		vec3_mult_vec(&_E2_2,&e2);
		vec3_mult(&_E2_2,2.0);
		vec3_t _e1_sq;
		vec3_copy(&_e1_sq,&e1);
		vec3_mult_vec(&_e1_sq,&e1);
		vec3_add_vec(&_E2_2,&_e1_sq);
		vec3_copy(&_E2_3,&e0);
		vec3_mult_vec(&_E2_3,&e1);
		vec3_mult(&_E2_3,2.0);
		vec3_copy(&_E2_4,&e0);
		vec3_mult_vec(&_E2_4,&e0);
		E_2[0] += vec3_sum(&_E2_0);
		E_2[1] += vec3_sum(&_E2_1);
		E_2[2] += vec3_sum(&_E2_2);
		E_2[3] += vec3_sum(&_E2_3);
		E_2[4] += vec3_sum(&_E2_4);
	}

	//scalar_array _a;
	//_a.resize(5);
	double _a[5];
	
	_a[4] = -E_2[1];
	_a[3] = (real_t)(4)*E_2[0] - (real_t)(2)*E_2[2];
	_a[2] = -(real_t)(3)*E_2[3] + (real_t)(3)*E_2[1];
	_a[1] = -(real_t)(4)*E_2[4] + (real_t)(2)*E_2[2];
	_a[0] = E_2[3];

	double at_sol[5];
	
	int num_sol = solve_polynomial(&*at_sol, &*_a, 5);
	double e[num_sol];
	scalar_array_clear(&*e, num_sol);
	
	double at[num_sol];
	
	if(COMPLICATED_ERROR_CMP)
	{
	  // get the error in a complicate way
	  scalar_array_clear(&*e, num_sol);
	  scalar_array_add(&*e, _a[0], num_sol);
	  
	  //at.clear();
	  //at.assign(at_sol.begin(),at_sol.end());
	  memcpy(&*at, &*at_sol, num_sol*sizeof(double));
	  
	  scalar_array_mult(&*at, _a[1], num_sol);
	  scalar_array_add_vec(&*e, &*at, num_sol);
	  
	  for(j=2; j<=4; j++)
	    {
	      //at.clear();
	      //at.assign(at_sol.begin(),at_sol.end());
		  memcpy(&*at, &*at_sol, num_sol*sizeof(double));
		  
	      scalar_array_pow(&*at, (real_t)(j), num_sol);
	      scalar_array_mult(&*at, _a[j], num_sol);
	      scalar_array_add_vec(&*e, &*at, num_sol);
	    }
	}
	else
	{
	  // Or in a fast one 
	  scalar_array_add(&*e, _a[4], num_sol);
	  for(j=3;j>0;j--)
	  {
	    // multiply with at & add a_ 
	    for(k=0;k<num_sol;k++)
		{
	      e[k] =  e[k]*at_sol[k] + _a[j] ;
		}
	  }
	}
	
	memcpy(&*at, &*at_sol, num_sol*sizeof(double));

	// get the angle al
	//scalar_array sa(at.begin(),at.end());
	double sa[num_sol];
	memcpy(&*sa, &*at, num_sol*sizeof(double));
	
	scalar_array_mult(&*sa, 2.0, num_sol);
	
	//scalar_array _ca1(at.begin(),at.end());
	double _ca1[num_sol];
	memcpy(&*_ca1, &*at, num_sol*sizeof(double));
	
	scalar_array_pow(&*_ca1,2.0, num_sol);
	scalar_array_add(&*_ca1,1.0, num_sol);
	
	//scalar_array ca(at.begin(),at.end());
	double ca[num_sol];
	memcpy(&*ca, &*at, num_sol*sizeof(double));
	
	scalar_array_pow(&*ca,2, num_sol);
	scalar_array_negate(&*ca, num_sol);
	scalar_array_add(&*ca,1.0, num_sol);
	scalar_array_div_vec(&*ca, &*_ca1, num_sol);
	scalar_array_div_vec(&*sa, &*_ca1, num_sol);
	
	double al[num_sol];
	scalar_array_atan2(&*al, &*sa, &*ca, num_sol);
	
	// check the sign of the derivative 
	scalar_array_mult(&*al, (real_t)(180./CONST_PI), num_sol);
	
	double _c_tMaxMin[num_sol];
	//_c_tMaxMin.resize(at.size());
	scalar_array_clear(&*_c_tMaxMin, num_sol);
	scalar_array_add(&*_c_tMaxMin, _a[1], num_sol);
	double _at[num_sol];
	//_at.clear();
	//_at.assign(at.begin(),at.end());
	memcpy(&*_at, &*at, num_sol*sizeof(double));
	
	scalar_array_mult(&*_at, _a[2], num_sol);
	scalar_array_mult(&*_at, 2.0, num_sol);
	scalar_array_add_vec(&*_c_tMaxMin, &*_at, num_sol);

	for(j=3; j<=4; j++)
	{
		memcpy(&*_at, &*at, num_sol*sizeof(double));
		
		scalar_array_pow(&*_at, (real_t)(j)-(real_t)(1.0), num_sol);
		scalar_array_mult(&*_at, _a[j], num_sol);
		scalar_array_mult(&*_at, (real_t)(j), num_sol);
		scalar_array_add_vec(&*_c_tMaxMin, &*_at, num_sol);
	}

	double tMaxMin[num_sol];
	double al_[num_sol];
	int al_idx = 0, a;
	
	memcpy(&*tMaxMin, &*_c_tMaxMin, num_sol*sizeof(double));
	
	for(i=0; i<num_sol; i++)
	{
		if(tMaxMin[i] > 0) al_[al_idx++] = al[i];
	}
	
	for(a=0; a<al_idx; a++)
	{
		vec3_t rpy;
		vec3_assign(&rpy, (real_t)0, (real_t)(al_[a] * CONST_PI / (real_t)(180)), (real_t)(0));
		mat33_t R, Ry_;
		rpyMat(&Ry_, &rpy);
		mat33_mult_mat2(&R, Rz, &Ry_);
		vec3_t t_opt;
		vec3_clear(&t_opt);

		for(i=0; i<n; i++)
		{
			mat33_t _m1, _eye3;
			mat33_eye(&_eye3);
			mat33_copy(&_m1, &V[i]);
			mat33_sub(&_m1, &_eye3);
			vec3_t _v1, _v2;
			vec3_mult_mat(&_v1, &R, &p[i]);
			vec3_mult_mat(&_v2, &_m1, &_v1);
			vec3_add_vec(&t_opt, &_v2);
		}

		vec3_t t_opt_;
		vec3_mult_mat(&t_opt_, &G, &t_opt);
		tnew[a] = t_opt_;
	}
	
	memcpy(al_ret, &*al_, al_idx*sizeof(double));
	
	*res_n = al_idx;
}

// =====================================================================================

void getRfor2ndPose_V_Exact(pose_t *sol, const vec3_t *v, const vec3_t *P,
					        const mat33_t R, const vec3_t t, const real_t DB, const int n, int *r_n)
{

	mat33_t RzN;
	decomposeR(&RzN, &R);
	mat33_t R_;
	mat33_mult_mat2(&R_, &R, &RzN);
	mat33_t RzN_tr;
	mat33_transpose(&RzN_tr, RzN);
	//vec3_array P_;
	vec3_t P_[n];
	
	vec3_array_mult_vec2(&*P_, &RzN_tr, P, n);
	
	vec3_t ang_zyx;
	rpyAng_X(&ang_zyx, &R_);
	
	vec3_t rpy;
	mat33_t Ry,Rz;
	vec3_assign(&rpy, 0, ang_zyx.v[1], 0);
	rpyMat(&Ry, &rpy);
	vec3_assign(&rpy,0,0,ang_zyx.v[2]);
	rpyMat(&Rz,&rpy);
	//scalar_array bl;
	//vec3_array Tnew;
	//double *bl = (double*)malloc(4 * sizeof(double));
	//vec3_t *Tnew = (vec3_t*)malloc(4 * sizeof(vec3_t));
	double bl[10];
	vec3_t Tnew[10];
	
	int res_n = 0;
	
	getRotationY_wrtT(&*bl, &*Tnew, v, &*P_, &t, &DB, &Rz, n, &res_n);
	
	// Estimate the Error for all solutions !
	scalar_array_div(bl, 180.0/CONST_PI, res_n);
	
	mat33_t V[n];
	
	int i, j;
	for(i=0; i<n; i++)
	{
		vec3_mul_vec3trans(&V[i], &v[i], &v[i]);
		mat33_div(&V[i], vec3trans_mul_vec3(&v[i],&v[i]));
	}
	
	mat33_t _m1;
	mat33_t _m2;
	vec3_t _v1;
	vec3_t _v2;

	for(j=0; j<res_n; j++)
	{
		mat33_clear(&Ry);
		vec3_assign(&rpy, 0, bl[j], 0);
		rpyMat(&Ry, &rpy);
		mat33_mult_mat2(&_m1, &Rz, &Ry);
		mat33_mult_mat2(&sol[j].R, &_m1, &RzN_tr);
		vec3_copy(&sol[j].t, &Tnew[j]);
		real_t E = 0;
		for(i=0; i<n; i++)
		{
			mat33_eye(&_m2);
			mat33_sub(&_m2, &V[i]);
			vec3_mult_mat(&_v1, &sol[j].R, &P[i]);
			vec3_add_vec(&_v1, &sol[j].t);
			vec3_mult_mat(&_v2, &_m2, &_v1);
			vec3_mult_vec(&_v2, &_v2);
			E += vec3_sum(&_v2);
		}
		sol[j].E = E;
	}
	
	*r_n = res_n;
}

// =====================================================================================

void get2ndPose_Exact(pose_t *sol, const vec3_t *v, const vec3_t *P,
					  const mat33_t R, const vec3_t t, const real_t DB, const int n, int *sol_n)
{
	vec3_t cent, _v1;
	
	vec3_t _va1[n];
	normRv2(&*_va1, v, n);
	
	vec3_array_mean(&_v1, &*_va1, n);
	normRv(&cent, &_v1);
	mat33_t Rim;
	vec3_clear(&_v1);
	_v1.v[2] = 1.0;
	
	GetRotationbyVector(&Rim, &_v1, &cent);
	
	vec3_t v_[n];
	vec3_array_mult_vec2(&*v_, &Rim, v, n);
	
	normRv2(&*_va1, &*v_, n);
	vec3_array_mean(&_v1, &*_va1, n);
	normRv(&cent, &_v1);
	mat33_t R_;
	vec3_t  t_;
	mat33_mult_mat2(&R_, &Rim, &R);
	vec3_mult_mat(&t_, &Rim, &t);
	
	int res_n = 0;
	getRfor2ndPose_V_Exact(sol, &*v_, P, R_, t_, DB, n, &res_n);
	
	// Re normalize it 
	mat33_t Rim_tr;
	mat33_transpose(&Rim_tr, Rim);
	
	int i;
	vec3_t _t;
	mat33_t _R;
	for(i=0; i<res_n; i++)
	{
		vec3_mult_mat(&_t, &Rim_tr, &sol[i].t);
		mat33_mult_mat2(&_R, &Rim_tr, &sol[i].R);

		vec3_copy(&sol[i].t, &_t);
		mat33_copy(&sol[i].R, &_R);
	}
	
	*sol_n = res_n;
}

// =====================================================================================
// function returns two solutions to the pose problem (1st has the lower Error !!!!!!
// ====================================================================================
void robust_pose(real_t *err, mat33_t *R , vec3_t *t, 
		 real_t *errb, mat33_t *Rb, vec3_t *tb,
		 const vec3_t *_model, const vec3_t *_iprts, const options_t _options, const int n)
{
 
  
  
	mat33_t Rlu_;
	vec3_t tlu_;
	int it1_, i;
	real_t obj_err1_;
	real_t img_err1_;
	
	vec3_t model[n];
	vec3_t iprts[n];
	memcpy(&*model, _model, n * sizeof(vec3_t));
	memcpy(&*iprts, _iprts, n * sizeof(vec3_t));
	
	options_t options;
	memcpy(&options, &_options, sizeof(options_t));

	mat33_clear(&Rlu_);
	vec3_clear(&tlu_);
	it1_ = 0;
	obj_err1_ = 0;
	img_err1_ = 0;
	
	objpose(&Rlu_, &tlu_, &it1_, &obj_err1_, &img_err1_, CALC_IMAGE_ERROR, &*model, &*iprts, options, n);

	pose_t sol[10];
	int sol_n = 0;
	
	get2ndPose_Exact(&*sol, &*iprts, &*model, Rlu_, tlu_, 0, n, &sol_n);
	
  // refine all poses 
  for(i=0; i<sol_n; i++)
  {
    mat33_copy(&options.initR, &sol[i].R);
	
    objpose(&Rlu_, &tlu_, &it1_, &obj_err1_, &img_err1_, CALC_IMAGE_ERROR, &*model, &*iprts, options, n);
	
    mat33_copy(&sol[i].PoseLu_R, &Rlu_);
    vec3_copy(&sol[i].PoseLu_t, &tlu_);
    sol[i].obj_err = obj_err1_;
  }
  
  //  from all solutions get the best 
  int min_err_idx = (-1);
  real_t min_err = MAX_FLOAT;
  for(i=0; i<sol_n; i++)
    if(sol[i].obj_err < min_err)
    {
      min_err = sol[i].obj_err;
      min_err_idx = i;
    }
  
  // Copy the best solution !!!
	if(min_err_idx >= 0)
	{
		mat33_copy(R, &sol[min_err_idx].PoseLu_R);
		vec3_copy(t, &sol[min_err_idx].PoseLu_t);
		*err = sol[min_err_idx].obj_err;
	}
	else
	{
		mat33_clear(R);
		vec3_clear(t);
		*err = MAX_FLOAT;
	}
//  from all solutions get the 2nd Best one !!
  int min_err_idx2 = min_err_idx;
  min_err = MAX_FLOAT;
  
  for(i=0; i < sol_n; i++)
  {
    if(sol[i].obj_err < min_err  && 
       (min_err_idx2 != i) )
    {
      min_err = sol[i].obj_err;
      min_err_idx = i;
    }
  }
  // Copy the 2nd best solution !!!
  if(min_err_idx >= 0)
  {
    mat33_copy(Rb, &sol[min_err_idx].PoseLu_R);
    vec3_copy(tb, &sol[min_err_idx].PoseLu_t);
    *errb = sol[min_err_idx].obj_err;
  }
  else
  {
    mat33_clear(R);
    vec3_clear(t);
    *err = MAX_FLOAT;
  }
  
}

// =====================================================================================
void robust_pose2(real_t *err, mat33_t *R, vec3_t *t,
		 const vec3_t *_model, const vec3_t *_iprts,
		 const options_t _options, const int n)
{
  real_t  errb;
  vec3_t  tb;
  mat33_t Rb;    
  
  // just call the function above 
  robust_pose(err, R, t, &errb, &Rb, &tb, _model, _iprts, _options, n);
}

// ----------------------------------------

