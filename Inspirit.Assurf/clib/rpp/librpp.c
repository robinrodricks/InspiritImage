#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#include "inc/librpp.h"
#include "inc/rpp.h"
#include "inc/rpp_vecmat.h"

//using namespace rpp;

void robustPlanarPose(rpp_float *err,
								 double *R,
								 double *t,
								 const rpp_float cc[2],
								 const rpp_float fc[2],
								 double *model,
								 double *iprts,
								 const int model_iprts_size,
								 double *R_init,
								 const int estimate_R_init,
								 const rpp_float epsilon,
								 const rpp_float tolerance,
								 const unsigned int max_iterations)
{
	vec3_t _model[model_iprts_size];
	vec3_t _iprts[model_iprts_size];
	
	//_model.resize(model_iprts_size);
	//_iprts.resize(model_iprts_size);

	mat33_t K, K_inv;
	mat33_eye(&K);
	K.m[0] = (real_t)fc[0];
	K.m[4] = (real_t)fc[1];
	K.m[2] = (real_t)cc[0];
	K.m[5] = (real_t)cc[1];

	mat33_inv(&K_inv, &K);

	int i;
	for(i=0; i<model_iprts_size; i++)
	{
		vec3_t _v,_v2;
		vec3_assign(&_v, (real_t)model[i*3+0], (real_t)model[i*3+1], (real_t)model[i*3+2]);
		_model[i] = _v;
		vec3_assign(&_v, (real_t)iprts[i*3+0], (real_t)iprts[i*3+1], (real_t)iprts[i*3+2]);
		vec3_mult_mat(&_v2, &K_inv, &_v);
		_iprts[i] = _v2;
	}

	options_t options;
	options.max_iter = max_iterations;
	options.epsilon = (real_t)(epsilon == 0 ? DEFAULT_EPSILON : epsilon);
	options.tol =     (real_t)(tolerance == 0 ? DEFAULT_TOL : tolerance);
	if(estimate_R_init == 1)
	{
		mat33_set_all_zeros(&options.initR);
	}
	else
	{
		/*mat33_assign(&options.initR,
					(real_t)R_init[0][0], (real_t)R_init[0][1], (real_t)R_init[0][2],
					(real_t)R_init[1][0], (real_t)R_init[1][1], (real_t)R_init[1][2],
					(real_t)R_init[2][0], (real_t)R_init[2][1], (real_t)R_init[2][2]);*/
	}

	real_t _err;
	mat33_t _R;
	vec3_t _t;

	robust_pose2(&_err, &_R, &_t, &*_model, &*_iprts, options, model_iprts_size);

	/*for(j=0; j<3; j++)
	{		
		R[j*3+0] = (rpp_float)_R.m[j][0];
		R[j*3+1] = (rpp_float)_R.m[j][1];
		R[j*3+2] = (rpp_float)_R.m[j][2];
		
		t[j] = (rpp_float)_t.v[j];
	}*/
	
	memcpy(&*R, &*_R.m, 9*sizeof(rpp_float));
	memcpy(&*t, &*_t.v, 3*sizeof(rpp_float));
	
	*err = (rpp_float)_err;
}
/*
double estimatePoseFromCorrespondences(double *Correspondences, double *data, const unsigned int num, const double iw, const double ih)
{
	const int initial_estimate_with_arGetInitRot = 1; // only for testing
	
	double err = 1e+20;
	rpp_mat R, R_init;
	rpp_vec t;
	unsigned int i, j;
	double center[2] = {iw*(double)-0.5, ih * (double)-0.5};
	
	rpp_vec ppos2d[num];
	rpp_vec ppos3d[num];
	
	const double model_z =  0;
	const double iprts_z =  1;

	for(i = 0; i < num; i++)
	{
		ppos2d[i][0] = Correspondences[i*6 + 2];//data[0];
		ppos2d[i][1] = Correspondences[i*6 + 3];//data[1];
		ppos2d[i][2] = iprts_z;

		ppos3d[i][0] = center[0] + Correspondences[i*6 + 4];//width*(double)0.5;
		ppos3d[i][1] = center[1] + Correspondences[i*6 + 5];//height*(double)0.5;
		ppos3d[i][2] = model_z;
	}

	const double cc[2] = { data[8], data[9] }; //{arCamera->mat[0][2],arCamera->mat[1][2]};
	const double fc[2] = { data[10], data[11] }; //{arCamera->mat[0][0],arCamera->mat[1][1]};
	
	robustPlanarPose(err,R,t,cc,fc,ppos3d,ppos2d,num,R_init, initial_estimate_with_arGetInitRot,0,0,0);
	
	for(i=0; i<3; i++)
	{
		//conv[i][3] = (double)t[i];
		data[i*4 + 3] = (double)t[i];
		for(j=0; j<3; j++)
		{
			//conv[i][j] = (double)R[i][j];
			data[i*4 + j] = (double)R[i][j];
		}
	}

	if(err > 1e+10) return(-1);
	
	return((double)err);
}

double estimatePlanarPose(double *data, double width, double height)
{
	const int initial_estimate_with_arGetInitRot = 1; // only for testing
	
	double err = 1e+20;
	rpp_mat R, R_init;
	rpp_vec t;
	int i, j;
	double center[2] = {0.0, 0.0};
	
	rpp_vec ppos2d[4];
	rpp_vec ppos3d[4];
	const unsigned int n_pts = 4;
	const double model_z =  0;
	const double iprts_z =  1;

	ppos2d[0][0] = data[0];
    ppos2d[0][1] = data[1];
	ppos2d[0][2] = iprts_z;
    ppos2d[1][0] = data[2];
    ppos2d[1][1] = data[3];
	ppos2d[1][2] = iprts_z;
    ppos2d[2][0] = data[4];
    ppos2d[2][1] = data[5];
	ppos2d[2][2] = iprts_z;
    ppos2d[3][0] = data[6];
    ppos2d[3][1] = data[7];
	ppos2d[3][2] = iprts_z;

	ppos3d[0][0] = center[0] - width*(double)0.5;
    ppos3d[0][1] = center[1] + height*(double)0.5;
	ppos3d[0][2] = model_z;
    ppos3d[1][0] = center[0] + width*(double)0.5;
    ppos3d[1][1] = center[1] + height*(double)0.5;
	ppos3d[1][2] = model_z;
    ppos3d[2][0] = center[0] + width*(double)0.5;
    ppos3d[2][1] = center[1] - height*(double)0.5;
	ppos3d[2][2] = model_z;
    ppos3d[3][0] = center[0] - width*(double)0.5;
    ppos3d[3][1] = center[1] - height*(double)0.5;
	ppos3d[3][2] = model_z;

	const double cc[2] = { data[8], data[9] }; //{arCamera->mat[0][2],arCamera->mat[1][2]};
	const double fc[2] = { data[10], data[11] }; //{arCamera->mat[0][0],arCamera->mat[1][1]};
	
	robustPlanarPose(err,R,t,cc,fc,ppos3d,ppos2d,n_pts,R_init, initial_estimate_with_arGetInitRot,0,0,0);
	
	for(i=0; i<3; i++)
	{
		//conv[i][3] = (double)t[i];
		data[i*4 + 3] = (double)t[i];
		for(j=0; j<3; j++)
		{
			//conv[i][j] = (double)R[i][j];
			data[i*4 + j] = (double)R[i][j];
		}
	}

	if(err > 1e+10) return(-1);
	
	return((double)err);
}*/