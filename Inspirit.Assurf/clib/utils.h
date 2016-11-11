#ifndef __MATH_UTILS__
#define __MATH_UTILS__

#include <math.h>

#include "kdtree/kdtree.h"

static const int IMG_BORDER = 10;
static const double NCC_THRESH = 0.80;

static const double pi = 3.14159;
static const double two_pi = 6.28318;
static const double pi_on_three = 1.04719667;

static const double INLIER_THRESHOLD_SQ = 100.0;
static const double PROBABILITY_REQUIRED = 0.99;
static const double SQRT2 = 1.4142135623730951;

// ----------------------------------------------

// Build fast ring orientation and index.

static const double fast_ring_x[16] = { 0,0.31622776601683794,0.7071067811865475,0.9486832980505138,1.0,0.9486832980505138,0.7071067811865475,0.31622776601683794,0,-0.31622776601683794,-0.7071067811865475,-0.9486832980505138,-1.0,-0.9486832980505138,-0.7071067811865475,-0.31622776601683794 };
static const double fast_ring_y[16] = { -1.0,-0.9486832980505138,-0.7071067811865475,-0.31622776601683794,0,0.31622776601683794,0.7071067811865475,0.9486832980505138,1.0,0.9486832980505138,0.7071067811865475,0.31622776601683794,0,-0.31622776601683794,-0.7071067811865475,-0.9486832980505138 };

static const int indX[16] = {3,3,2,1,0,-1,-2,-3,-3,-3,-2,-1,0,1,2,3};
static const int indY[16] = {0,1,2,3,3,3,2,1,0,-1,-2,-3,-3,-3,-2,-1};

// ----------------------------------------------

typedef struct _IPoint{
	int index, localIndex, refIndex, pos;
	int x, y, dx, dy, pyrLevel, matched;
	double score;
	double orientation;
	int sum, sqsum, mean;
	double stdev;
	double *descriptor;
	int sampled;
	unsigned char sample[256];
}IPoint;

typedef struct _IPointMatch{
	IPoint *first;
	IPoint *second;
	IPoint *prev;
	double confidence;
	double normConfidence;
	int lk_tracked, tracked;
	int wasGood;
}IPointMatch;

typedef struct _RefObject{
	int index, pointsCount, matchedPointsCount, prevMatchedPointsCount;
	int width, height;
	double poseError;
	IPoint *points;
	IPointMatch *matches;
	double *descriptors;
	double homography[9];
	double pose[12];
	VlKDForest *kdf;
}RefObject;

int compareIPoint(const void *vp1, const void *vp2);
int compareIPointMatch(const void *vp1, const void *vp2);
void perpendicularRegression(IPointMatch *matches, const int count, IPointMatch *resultMatches);
void perpendicularRegressionIdx(IPointMatch *matches, const int count, int *resultMatches);
void getHomographyPoints(IPointMatch *matches, const int count, double *scrPts, double *refPts);
void getHomographyPoints2(IPointMatch *matches, const int count, double *scrPts, double *refPts);
void sortMatchesByObjects(RefObject *objMap, IPointMatch *matches, const int objNum, const int matchNum);
double bilinear_interpolation(const unsigned char *arr, const int w, const double x, const double y);

// ----------------------------------------------

int filterOutliersByAngle(IPointMatch *matches, const int matchedCount);
int filterOutliersByLines(IPointMatch *matches, const int matchedCount);

// ----------------------------------------------

inline int IMIN(register int a, register int b);
inline int IMAX(register int a, register int b);
inline double FMAX(register double a, register double b);
inline double IGAUSSIAN(register int x, register int y, register double sig);
inline double FGAUSSIAN(register double x, register double y, register double sig);
inline double ANGLE(const double X, const double Y);
inline double fast_atan2(const double y, const double x);
inline int dRound(register double dbl);
inline double dSquare(register double dbl);
inline void sin_cos(const double a, double *s, double *c);
inline double fast_sqrt(const double x);
inline int iabs(register int value);
inline int iSquare(register int dbl);
inline int testSideOfLine(const int Ax, const int Ay, const int Bx, const int By, const int Cx, const int Cy);
inline double getCoterminalAngle(double angle);

#endif
