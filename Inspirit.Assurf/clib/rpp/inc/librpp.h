
#include "rpp_types.h"

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
								 

//double estimatePlanarPose(double *data, double width, double height);
//extern "C" double estimatePoseFromCorrespondences(double *Correspondences, double *data, const unsigned int num, const double iw, const double ih);

/*

	[OUTPUT]

	err: squared reprojection error
	R:   rotation matrix (iprts[n] = R*model[n]+t)
	t:   translation vector (iprts[n] = R*model[n]+t)

	[INPUT]
    
	cc:    camera's principal point [x,y]
	fc:    camera's focal length    [x,y]
	model: 3d points [x,y,z]
	iprts: 2d projections [x,y,1]
    
    model_iprts_size:   number of 2d/3d point correspondences
    R_init:             initial estimate of the rotation matrix R
    estimate_R_init:    when true, the estimate in R_init is ignored
    epsilon*:           see below (default: 1E-8)
    tolerance*:         see below (default: 1E-5)
	max_iterations*:    max. number of iterations (0 = infinite)
    
    *) the following code fragment illustrates the use of epsilon,
       tolerance and max_iterations:

    while((ABS(( old_err - new_err ) / old_err) > tolerance) && 
          ( new_err > epsilon ) &&
		  ( max_iterations == 0 || iterations < max_iterations ))
    {
        NEW ITERATION
    }
          
*/
