#include <math.h>
#include "utils.h"

// ----------------------------------------------

typedef double rpp_float;

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
								 const unsigned int max_iterations);

// ----------------------------------------------

double estimatePoseFromMatches(double *camera_info, IPointMatch *matches, double model_matrix[12], 
								const int matchesCount, const int width, const int height)
{
	const int initial_estimate_with_arGetInitRot = 1; // only for testing
	
	double err = 1e+20;
	double R[3*3];
	double R_init[3*3];
	double t[3];
	int i, j;
	double center[2] = {(double)width * -0.5, (double)height * -0.5};
	
	double ppos2d[matchesCount*3];
	double ppos3d[matchesCount*3];
	
	const double model_z =  0.0;
	const double iprts_z =  1.0;

	for(i = 0; i < matchesCount; ++i)
	{
		ppos2d[i*3+0] = matches[i].first->x;
		ppos2d[i*3+1] = matches[i].first->y;
		ppos2d[i*3+2] = iprts_z;

		ppos3d[i*3+0] = center[0] + matches[i].second->x;
		ppos3d[i*3+1] = center[1] + matches[i].second->y;
		ppos3d[i*3+2] = model_z;
	}

	const double cc[2] = { camera_info[0], camera_info[1] };
	const double fc[2] = { camera_info[2], camera_info[3] };
	
	robustPlanarPose(&err, &*R, &*t, cc, fc, &*ppos3d, &*ppos2d, matchesCount, &*R_init, initial_estimate_with_arGetInitRot,0,0,1000);
	
	for(i=0; i<3; i++)
	{
		model_matrix[i*4 + 3] = (double)t[i];
		for(j=0; j<3; j++)
		{
			model_matrix[i*4 + j] = (double)R[i*3+j];
		}
	}

	if(err > 1e+10) return 1e+10;
	
	return((double)err);
}