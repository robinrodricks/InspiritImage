#ifndef MYPYRLK_H
#define MYPYRLK_H

#include "utils.h"

#define DBL_EPSILON	2.2204460492503131E-16
#define CV_LKFLOW_INITIAL_GUESSES   4
#define CV_LKFLOW_GET_MIN_EIGENVALS 8
#define CV_TERMCRIT_ITER    1
#define CV_TERMCRIT_NUMBER  CV_TERMCRIT_ITER
#define CV_TERMCRIT_EPS     2

#define uchar unsigned char

typedef struct _CvPoint2D32f{
	double x, y;
}CvPoint2D32f;

typedef struct _CvSize{
	int width, height;
}CvSize;

typedef struct _CvPoint{
	int x, y;
}CvPoint;

typedef struct _CvRect{
    int x;
    int y;
    int width;
    int height;
}CvRect;

typedef struct _CvTermCriteria{
    int    type;
    int    max_iter;
    double epsilon;
}CvTermCriteria;

void initOpticalFlow( uchar* arrA[3], uchar* arrB[3], 
						const int sizeW[3], const int sizeH[3], int patch_size, int flags);

void myCalcOpticalFlowPyrLK( uchar* arrA[3], uchar* arrB[3],
						const int sizeW[3], const int sizeH[3],
                        const CvPoint2D32f * featuresA, CvPoint2D32f * featuresB,
                        int count, int patch_size, char *status, double *error,
                        int flags );


#endif
