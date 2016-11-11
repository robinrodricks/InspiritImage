#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "utils.h"


int compareIPoint(const void *vp1, const void *vp2)
{
     const IPoint *pkt_a = (const IPoint *)vp1;
     const IPoint *pkt_b = (const IPoint *)vp2;

     if (pkt_a->score <  pkt_b->score) return 1;
     if (pkt_a->score >  pkt_b->score) return -1;
     return 0;
}

int compareIPointMatch(const void *vp1, const void *vp2)
{
     const IPointMatch *pkt_a = (const IPointMatch *)vp1;
     const IPointMatch *pkt_b = (const IPointMatch *)vp2;

     if (pkt_a->confidence <  pkt_b->confidence) return -1;
     if (pkt_a->confidence >  pkt_b->confidence) return 1;
     return 0;
}

// ----------------------------------------------

inline int IMIN(register int a, register int b)
{
	return (((a) < (b)) ? (a) : (b));
}

inline int IMAX(register int a, register int b)
{
	return (((a) < (b)) ? (b) : (a));
}

inline double FMAX(register double a, register double b)
{
	return (((a) < (b)) ? (b) : (a));
}

inline double IGAUSSIAN(register int x, register int y, register double sig)
{
	return (((double)1.0 / (two_pi*sig*sig) * exp( -(x*x+y*y) / ((double)2.0*sig*sig))));
}
inline double FGAUSSIAN(register double x, register double y, register double sig)
{
	return (((double)1.0 / (two_pi*sig*sig) * exp( -(x*x+y*y) / ((double)2.0*sig*sig))));
}

inline double ANGLE(const double X, const double Y)
{
	if(X > 0 && Y >= 0) return atan(Y/X);

	if(X < 0 && Y >= 0) return pi - atan(-Y/X);

	if(X < 0 && Y < 0) return pi + atan(Y/X);

	if(X > 0 && Y < 0) return two_pi - atan(-Y/X);

	return 0.0;
}

inline double fast_atan2(const double y, const double x)
{
	double angle, r;
	const double sign = y < 0 ? -1.0 : 1.0;
	double abs_y = y*sign + 2.220446049250313e-16;

	if (x >= 0) {
		r = (x - abs_y) / (x + abs_y);
		angle = (double)0.785398164;
	} else {
		r = (x + abs_y) / (abs_y - x);
		angle = (double)2.356194491;
	}
	angle += ((double)0.1821 * r * r - (double)0.9675) * r;
	return angle * sign;
}

inline int iabs(register int value)
{
	return ((value < 0) ? (-value) : (value));
}

inline int dRound(register double dbl)
{
	return (int) (dbl+0.5);
}

inline double dSquare(register double dbl)
{
	return (dbl * dbl);
}
inline int iSquare(register int dbl)
{
	return (dbl * dbl);
}

inline void sin_cos(const double a, double *s, double *c)
{
	double angleTmp = a * 0.6366197723675814;

	int qi = (angleTmp + (double)0.5 - (double)(angleTmp < 0 ? 1 : 0));
	angleTmp -= qi;

	double xq = angleTmp * angleTmp;

	//more precision
	double f0 = (double)1.0 + xq * (xq * (0.25360671639164339 - 0.020427240364907607 * xq) - 1.2336979844380824);
	double f1 = angleTmp * (1.5707963267948966 + xq * (  xq * (0.079679708649230657 - 0.0046002309092153379 * xq) - 0.64596348437163809));

	double qi1 = qi & 1;
	double qi2 = qi & 2;
	double qi1_2 = (double)1.0 - qi2;
	qi2 = qi1 * qi1_2;
	qi1 = ((double)1.0 - qi1) * qi1_2;

	*c = qi1 * f0 - qi2 * f1;
	*s = qi1 * f1 + qi2 * f0;
}

inline double fast_sqrt(const double x) 
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

// Return the coterminal angle between [0;2*PI].
inline double getCoterminalAngle(double angle)
{
	while (angle > two_pi) angle -= two_pi;
	while (angle < 0.0) angle += two_pi;
	return angle;
}

inline int testSideOfLine(const int Ax, const int Ay, const int Bx, const int By, const int Cx, const int Cy)
{
	// 0 - CxCy is on the line
	// +1 - CxCy is on the right
	// -1 - CxCy is on the left
	int s = (Bx - Ax) * (Cy - Ay) - (By - Ay) * (Cx - Ax);
	if(s < 0) return -1;
	if(s > 0) return +1;
	return 0;
}

int solve_deg2(double a, double b, double c, double *x1, double *x2)
{
	double delta = b * b - (double)4.0 * a * c;

	if (delta < 0) return 0;

	double inv_2a = (double)0.5 / a;

	if (delta == 0)
	{
		*x1 = -b * inv_2a;
		*x2 = *x1;
		return 1;
	}

	double sqrt_delta = fast_sqrt(delta);
	*x1 = (-b + sqrt_delta) * inv_2a;
	*x2 = (-b - sqrt_delta) * inv_2a;
	return 2;
}

int filterOutliersByAngle(IPointMatch *matches, const int matchedCount)
{
	int i, j, k;
	const int nbins = 360 / 4;
	double histogram[nbins];
	int rotdiff[matchedCount];
	const double div = 57.295827908797776;
	int r0, r1, ar;
	
	for(i=0; i<nbins;i++)histogram[i]=0;
	
	for(i = 0; i < matchedCount; i++)
	{
		const IPointMatch *mtch = &matches[i];
		
		r0 = mtch->first->orientation * div;
		r1 = mtch->second->orientation * div;
		ar = r1 - r0;
		if(ar < 0) ar += 360;
		histogram[ (ar)/4 ] += (double)1.1-mtch->normConfidence;
		rotdiff[i] = ar;
	}
	
	int index = 0;
	double maxBin = histogram[index];
    for (i = 1; i < nbins; i++) 
	{
		if (maxBin < histogram[i]) 
		{
			index = i;
			maxBin = histogram[i];
		}
    }
	
	index *= 4;
	
	k = 0;
	j = 0;
	const int thresh = 10*10;
	for(i = 0; i < matchedCount; i++)
	{
		if( iSquare(rotdiff[i] - index) <= thresh )
		{
			matches[k++] = matches[ i ];
		}
	}
	return k;
}

int filterOutliersByLines(IPointMatch *matches, const int matchedCount)
{
	if(matchedCount > 15)
	{
		qsort(matches, matchedCount, sizeof(IPointMatch), compareIPointMatch);
		
		int i, j;
		int index = IMIN(matchedCount >> 1, 30);
		int half = matchedCount >> 1;
		int imap[matchedCount];
		
		for(i=0;i<matchedCount;i++)imap[i]=0;
		
		for(i = 0; i < index; i++)
		{
			const int ind1 = i << 1;
			const int ind2 = ind1 + 1;
			const IPointMatch *m0 = &matches[ind1];
			const IPointMatch *m1 = &matches[ind2];
			int fail = 0;
			for(j = 0; j < matchedCount; j++)
			{
				if( ind1 == j || ind2 == j ) continue;
				const IPointMatch *m2 = &matches[j];
				if( testSideOfLine(m0->first->x, m0->first->y, m1->first->x, m1->first->y, m2->first->x, m2->first->y) != 
					testSideOfLine(m0->second->x, m0->second->y, m1->second->x, m1->second->y, m2->second->x, m2->second->y) )
				{
					if(++fail > half-2)
					{
						imap[ind1] += 1;
						imap[ind2] += 1;
						break;
					}
				}
			}
		}
		index = 0;
		for(i = 0; i < matchedCount; i++)
		{
			if(imap[i] < 1)
			{
				matches[index++] = matches[i];
			}
		}
		return index;
	} 
	else 
	{
		return matchedCount;
	}
}

void get2ClosestPoints(const int px, const int py, IPointMatch *matches, const int count, IPointMatch *result, int *matchMask)
{
	int i;
	int d1 = 2147483647, d2 = 2147483647;
	int ind1 = 0, ind2 = 0;
	
	for(i = 0; i < count; i++)
	{
		if(matchMask[i] == 1) continue;
		
		const int dist = iSquare(px - matches[i].first->x) + iSquare(py - matches[i].first->y);
		if(dist < d1)
		{
			d2 = d1;
			d1 = dist;
			ind2 = ind1;
			ind1 = i;
		}
		else if(dist < d2)
		{
			d2 = dist;
			ind2 = i;
		}
	}
	result[0] = matches[ind1];
	result[1] = matches[ind2];
	matchMask[ind1] = 1;
	matchMask[ind2] = 1;
}

void get2ClosestPointsIdx(const int px, const int py, IPointMatch *matches, const int count, int *result, int *matchMask)
{
	int i;
	int d1 = 2147483647, d2 = 2147483647;
	int ind1 = 0, ind2 = 0;
	
	for(i = 0; i < count; i++)
	{
		if(matchMask[i] == 1) continue;
		
		const int dist = iSquare(px - matches[i].first->x) + iSquare(py - matches[i].first->y);
		if(dist < d1)
		{
			d2 = d1;
			d1 = dist;
			ind2 = ind1;
			ind1 = i;
		}
		else if(dist < d2)
		{
			d2 = dist;
			ind2 = i;
		}
	}
	result[0] = ind1;
	result[1] = ind2;
	matchMask[ind1] = 1;
	matchMask[ind2] = 1;
}

void getHomographyPoints(IPointMatch *matches, const int count, double *scrPts, double *refPts)
{
	int i;
	int sx = 0, sy = 0, sxy = 0, sxx = 0;
	int minx = 1000, maxx = 0;
	int miny = 1000, maxy = 0;
	
	int result_msk[count];
	
	for(i = 0; i < count; i++)
	{
		const IPoint *pt = matches[i].first;
		sx += pt->x;
		sy += pt->y;
		sxy += pt->x * pt->y;
		sxx += iSquare(pt->x);
		minx = IMIN(minx, pt->x);
		maxx = IMAX(maxx, pt->x);
		miny = IMIN(miny, pt->y);
		maxy = IMAX(maxy, pt->y);
		
		result_msk[i] = 0;
	}
	
	double b = (double)( count * sxy - (sx * sy) ) / (double)( count * sxx - sx*sx );
	double a = (double)((double)sy - b * (double)sx) / (double)count;
	double ib = -1.0 / b;
	// y = a + b * x;
	// y/b - a/b = x;
	// perp
	// y = c + (-1/b) * x;
	// c = y - ((-1/b) * x)
	
	int lx = minx;
	int ly = a + b * lx;
	int rx = maxx;
	int ry = a + b * rx;
	
	int cx = (lx + rx) * 0.5;
	int cy = (ly + ry) * 0.5;
	
	double c = (double)cy - (ib * (double)cx);
	
	int ttx = (double)miny / ib - c / ib;
	int tty = miny;

	int btx = (double)maxy / ib - c / ib;
	int bty = maxy;
	
	IPointMatch result_pt[8];
	
	get2ClosestPoints(ttx, tty, matches, count, &*result_pt, &*result_msk); // a - b // 0 - 1
	get2ClosestPoints(btx, bty, matches, count, &*(result_pt+2), &*result_msk); // c - d // 2 - 3
	get2ClosestPoints(lx, ly, matches, count, &*(result_pt+4), &*result_msk); // e - f // 4 - 5
	get2ClosestPoints(rx, ry, matches, count, &*(result_pt+6), &*result_msk); // g - h // 6 - 7
	
	const int idx[16][4] = {
		{0, 6, 2, 4},
		{0, 7, 2, 4},
		{0, 6, 3, 4},
		{0, 7, 3, 4},
		{0, 7, 3, 5},
		{0, 6, 2, 5},
		{0, 6, 3, 5},
		{0, 7, 2, 5},
		
		{1, 6, 2, 4},
		{1, 7, 2, 4},
		{1, 6, 3, 4},
		{1, 7, 3, 4},
		{1, 7, 3, 5},
		{1, 6, 2, 5},
		{1, 6, 3, 5},
		{1, 7, 2, 5}
	};
	
	int j, k = 0, l = 0;
	for(i = 0; i < 16; i++)
	{
		for(j = 0; j < 4; j++)
		{
			const IPoint *pt0 = result_pt[ idx[i][j] ].first;
			const IPoint *pt1 = result_pt[ idx[i][j] ].second;
			scrPts[k++] = pt0->x;
			scrPts[k++] = pt0->y;
			
			refPts[l++] = pt1->x;
			refPts[l++] = pt1->y;
		}
	}
}

void perpendicularRegression(IPointMatch *matches, const int count, IPointMatch *resultMatches)
{
	int i;
	double meanx = 0.0, meany = 0.0;
	int minx = 1000, maxx = 0;
	int miny = 1000, maxy = 0;
	
	int result_msk[count];

	for(i = 0; i < count; i++)
	{
		const IPoint *pt = matches[i].first;
		meanx += pt->x;
		meany += pt->y;

		minx = IMIN(minx, pt->x);
		maxx = IMAX(maxx, pt->x);
		miny = IMIN(miny, pt->y);
		maxy = IMAX(maxy, pt->y);
		
		result_msk[i] = 0;
	}
	meanx /= (double)count;
	meany /= (double)count;

	double sxy = 0.0, sx2y2 = 0.0;
	for(i = 0; i < count; i++)
	{
		const IPoint *pt = matches[i].first;
		double nx = pt->x - meanx;
		double ny = pt->y - meany;
		sxy += nx * ny;
		sx2y2 += nx*nx - ny*ny;
	}

	double b = sx2y2 / sxy;

	// tan(q)^2 + b*tan(q) - 1 = 0
	// a*x^2 + bx + c = 0;
	// x = ( -b + Math.sqrt(b*b - 4*a*c) ) / 2*a
	// x = ( -b - Math.sqrt(b*b - 4*a*c) ) / 2*a
	const double toDeg = 57.29577951308232;

	double a_1 = 0.0;// = ( -b + sqrt(b * b - 4*-1) ) / 2;
	double a_2 = 0.0;// = ( -b - sqrt(b * b - 4*-1) ) / 2;
	
	solve_deg2(1.0, b, -1.0, &a_1, &a_2);
	
	double a1 = atan( a_1 ) * toDeg;
	double a2 = atan( a_2 ) * toDeg;
	if(a1 < 0) a1+=360;
	if(a2 < 0) a2+=360;
	if(a1 > 360) a1-=360;
	if(a2 > 360) a2-=360;

	double c1 = meany - (a_1 * meanx);
	double c2 = meany - (a_2 * meanx);
	
	int ttx, tty, btx, bty, lx, ly, rx, ry;

	if((a2 >= 45 && a2 <= 90+45) || (a2 >= 180+45 && a2 <= 270+45))
	{
		ttx = ((double)miny - c2) / a_2;
		tty = miny;

		btx = ((double)maxy - c2) / a_2;
		bty = maxy;
	} else {
		ttx = minx;
		tty = c2 + a_2 * minx;

		btx = maxx;
		bty = c2 + a_2 * maxx;
	}

	if((a1 >= 45 && a1 <= 90+45) || (a1 >= 180+45 && a1 <= 270+45))
	{
		lx = ((double)miny - c1) / a_1;
		ly = miny;

		rx = ((double)maxy - c1) / a_1;
		ry = maxy;
	} else {
		lx = minx;
		ly = c1 + a_1 * minx;

		rx = maxx;
		ry = c1 + a_1 * maxx;
	}
	
	get2ClosestPoints(ttx, tty, matches, count, &*resultMatches, &*result_msk); // a - b // 0 - 1
	get2ClosestPoints(btx, bty, matches, count, &*(resultMatches+2), &*result_msk); // c - d // 2 - 3
	get2ClosestPoints(lx, ly, matches, count, &*(resultMatches+4), &*result_msk); // e - f // 4 - 5
	get2ClosestPoints(rx, ry, matches, count, &*(resultMatches+6), &*result_msk); // g - h // 6 - 7
}

void perpendicularRegressionIdx(IPointMatch *matches, const int count, int *resultMatches)
{
	int i;
	double meanx = 0.0, meany = 0.0;
	int minx = 1000, maxx = 0;
	int miny = 1000, maxy = 0;
	
	int result_msk[count];

	for(i = 0; i < count; i++)
	{
		const IPoint *pt = matches[i].first;
		meanx += pt->x;
		meany += pt->y;

		minx = IMIN(minx, pt->x);
		maxx = IMAX(maxx, pt->x);
		miny = IMIN(miny, pt->y);
		maxy = IMAX(maxy, pt->y);
		
		result_msk[i] = 0;
	}
	meanx /= (double)count;
	meany /= (double)count;

	double sxy = 0.0, sx2y2 = 0.0;
	for(i = 0; i < count; i++)
	{
		const IPoint *pt = matches[i].first;
		double nx = pt->x - meanx;
		double ny = pt->y - meany;
		sxy += nx * ny;
		sx2y2 += nx*nx - ny*ny;
	}

	double b = sx2y2 / sxy;

	// tan(q)^2 + b*tan(q) - 1 = 0
	// a*x^2 + bx + c = 0;
	// x = ( -b + Math.sqrt(b*b - 4*a*c) ) / 2*a
	// x = ( -b - Math.sqrt(b*b - 4*a*c) ) / 2*a
	const double toDeg = 57.29577951308232;

	double a_1 = 0.0;// = ( -b + sqrt(b * b - 4*-1) ) / 2;
	double a_2 = 0.0;// = ( -b - sqrt(b * b - 4*-1) ) / 2;
	
	solve_deg2(1.0, b, -1.0, &a_1, &a_2);
	
	double a1 = atan( a_1 ) * toDeg;
	double a2 = atan( a_2 ) * toDeg;
	if(a1 < 0) a1+=360;
	if(a2 < 0) a2+=360;
	if(a1 > 360) a1-=360;
	if(a2 > 360) a2-=360;

	double c1 = meany - (a_1 * meanx);
	double c2 = meany - (a_2 * meanx);
	
	int ttx, tty, btx, bty, lx, ly, rx, ry;

	if((a2 >= 45 && a2 <= 90+45) || (a2 >= 180+45 && a2 <= 270+45))
	{
		ttx = ((double)miny - c2) / a_2;
		tty = miny;

		btx = ((double)maxy - c2) / a_2;
		bty = maxy;
	} else {
		ttx = minx;
		tty = c2 + a_2 * minx;

		btx = maxx;
		bty = c2 + a_2 * maxx;
	}

	if((a1 >= 45 && a1 <= 90+45) || (a1 >= 180+45 && a1 <= 270+45))
	{
		lx = ((double)miny - c1) / a_1;
		ly = miny;

		rx = ((double)maxy - c1) / a_1;
		ry = maxy;
	} else {
		lx = minx;
		ly = c1 + a_1 * minx;

		rx = maxx;
		ry = c1 + a_1 * maxx;
	}
	
	get2ClosestPointsIdx(ttx, tty, matches, count, &*resultMatches, &*result_msk); // a - b // 0 - 1
	get2ClosestPointsIdx(btx, bty, matches, count, &*(resultMatches+2), &*result_msk); // c - d // 2 - 3
	get2ClosestPointsIdx(lx, ly, matches, count, &*(resultMatches+4), &*result_msk); // e - f // 4 - 5
	get2ClosestPointsIdx(rx, ry, matches, count, &*(resultMatches+6), &*result_msk); // g - h // 6 - 7
}

void getHomographyPoints2(IPointMatch *matches, const int count, double *scrPts, double *refPts)
{
	IPointMatch result_pt[8];
	perpendicularRegression(matches, count, &*result_pt);
	
	const int idx[16][4] = {
		{0, 6, 2, 4},
		{0, 7, 2, 4},
		{0, 6, 3, 4},
		{0, 7, 3, 4},
		{0, 7, 3, 5},
		{0, 6, 2, 5},
		{0, 6, 3, 5},
		{0, 7, 2, 5},
		
		{1, 6, 2, 4},
		{1, 7, 2, 4},
		{1, 6, 3, 4},
		{1, 7, 3, 4},
		{1, 7, 3, 5},
		{1, 6, 2, 5},
		{1, 6, 3, 5},
		{1, 7, 2, 5}
	};
	
	int i, j, k = 0, l = 0;
	for(i = 0; i < 16; i++)
	{
		for(j = 0; j < 4; j++)
		{
			const IPoint *pt0 = result_pt[ idx[i][j] ].first;
			const IPoint *pt1 = result_pt[ idx[i][j] ].second;
			scrPts[k++] = pt0->x;
			scrPts[k++] = pt0->y;
			
			refPts[l++] = pt1->x;
			refPts[l++] = pt1->y;
		}
	}
}

double bilinear_interpolation(const unsigned char *arr, const int w, const double x, const double y)
{
	int mnx = (int)floor( x );
	int mny = (int)floor( y );
	int mxx = (int) ceil( x );
	int mxy = (int) ceil( y );

	double alfa = (double)mxx - x;
	double beta = (double)mxy - y;

	if( alfa < 0.001 ) alfa = 0;
	if( beta < 0.001 ) beta = 0;

	int mnyw = mny * w;
	int mxyw = mxy * w;

	if( alfa < 0.001 ) return (double)(beta * arr[mnyw+mxx] + ((double)1.0-beta) * arr[mxyw+mxx]);
	if( alfa > 0.999 ) return (double)(beta * arr[mnyw+mnx] + ((double)1.0-beta) * arr[mxyw+mnx]);
	if( beta < 0.001 ) return (double)(alfa * arr[mxyw+mnx] + ((double)1.0-alfa) * arr[mxyw+mxx]);
	if( beta > 0.999 ) return (double)(alfa * arr[mnyw+mnx] + ((double)1.0-alfa) * arr[mnyw+mxx]);

	return (double)( beta*(alfa * arr[mnyw+mnx] + ((double)1.0-alfa)*arr[mnyw+mxx] )
	            +((double)1.0-beta)*(alfa * arr[mxyw+mnx] + ((double)1.0-alfa)*arr[mxyw+mxx] ) );
}

void sortMatchesByObjects(RefObject *objMap, IPointMatch *matches, const int objNum, const int matchNum)
{
	int i;
	RefObject *obj;
	
	for(i = 0; i < objNum; i++)
	{
		objMap[i].prevMatchedPointsCount = objMap[i].matchedPointsCount;
		objMap[i].matchedPointsCount = 0;
	}
	
	for(i = 0; i < matchNum; i++)
	{
		IPointMatch *m = &matches[i];
		obj = &objMap[ m->second->refIndex ];
		memcpy(obj->matches + obj->matchedPointsCount, m, sizeof(IPointMatch));
		obj->matchedPointsCount++;
	}
}