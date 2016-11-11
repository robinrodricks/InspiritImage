
#include "libfast/fast.h"

int shiTomasiWin[196];

void calcShiTomasiWin(const int stride, const int nHalfBoxSize);
double FindShiTomasiScoreAtPoint(const unsigned char *image, const int stride, const xy *irCenter);
void fastRotationEstimation(const unsigned char *image, IPoint *ip);
int supressNeighborPoints(IPoint *pointsData, const int cnt, const double maxDist);


int detectPointsLevel(const int pyrLevel, const int width, const int height, IPoint *pointsData,
						const int corn_thresh, const double tomasi_thresh, const int maxPoints, const int imgBorder, const int maxPointsPool)
{
	int cNum, i, cnt = 0;
	
	const unsigned char *image = pyr_img[pyrLevel];
	
	xy *corners;
	corners = fast9_detect_nonmax(image, width, height, width, corn_thresh, &cNum);
	
	if(cNum > maxPointsPool)
	{
		cNum = maxPointsPool;
	}
	
	calcShiTomasiWin(width, 3);
	
	for(i=0; i < cNum; i++)
	{
		const xy *corn = &corners[i];
		
		if( useMask == 1 && img_mask[(corn->y<<pyrLevel) * pyr_img1_width + (corn->x<<pyrLevel)] == 0 ) continue;
		
		if(corn->x > imgBorder && corn->y > imgBorder && corn->x < width - imgBorder && corn->y < height - imgBorder)
		{
			const double score = FindShiTomasiScoreAtPoint(image, width, corn);
			
			if(score > tomasi_thresh)
			{
				IPoint *ipt = &pointsData[cnt++];
				
				ipt->x = corn->x<<pyrLevel;
				ipt->y = corn->y<<pyrLevel;
				ipt->dx = ipt->dy = 0;
				ipt->pos = corn->y * width + corn->x;
				ipt->pyrLevel = pyrLevel;
				ipt->score = score;
				ipt->sampled = 0;
				ipt->matched = 0;
			}
		}
	}
	
	free(corners);
	
	if(cnt > maxPoints)
	{
		qsort(pointsData, cnt, sizeof(IPoint), compareIPoint);
		cnt = maxPoints;
	}
	
	for(i = 0; i < cnt; i++)
	{		
		fastRotationEstimation(image, &pointsData[i]);
	}
	
	return cnt;
}

int detectPointsLevel2(const int numLevels, IPoint *pointsData, const int corn_thresh, 
						const double tomasi_thresh, const int maxPoints, const int imgBorder, const int maxPointsPool)
{
	int cNum, i, j, k, cnt = 0;
	xy *corners;
	IPoint *ipt;
	
	for(j = 0; j < numLevels; j++)
	{
		const unsigned char *image = pyr_img[j];
		const int width = pyr_img1_width >> j;
		const int height = pyr_img1_height >> j;
		k = 0;
		
		corners = fast9_detect_nonmax(image, width, height, width, corn_thresh, &cNum);
		
		if(cNum > maxPointsPool-cnt)
		{
			cNum = maxPointsPool-cnt;
		}
		
		calcShiTomasiWin(width, 3);
		
		for(i = 0; i < cNum; i++)
		{
			const xy *corn = &corners[i];
			
			if( useMask == 1 && img_mask[(corn->y<<j) * pyr_img1_width + (corn->x<<j)] == 0 ) continue;
			
			if(corn->x > imgBorder && corn->y > imgBorder && corn->x < width - imgBorder && corn->y < height - imgBorder)
			{
				const double score = FindShiTomasiScoreAtPoint(image, width, corn);
				
				if(score > tomasi_thresh)
				{
					ipt = &pointsData[cnt++];
					
					ipt->x = corn->x<<j;
					ipt->y = corn->y<<j;
					ipt->dx = ipt->dy = 0;
					ipt->pos = corn->y * width + corn->x;
					ipt->pyrLevel = j;
					ipt->score = score;
					ipt->sampled = 0;
					ipt->matched = 0;
					ipt->index = cnt;
					
					k++;
				}
			}
		}
		
		free(corners);
		
		//i = k;
		//k = supressNeighborPoints(pointsData+cnt-k, k, 9);
		//cnt -= i - k;
	}
	
	cnt = supressNeighborPoints(pointsData, cnt, 13.0);
	
	if(cnt > maxPoints)
	{
		//qsort(pointsData, cnt, sizeof(IPoint), compareIPoint);
		cnt = maxPoints;
	}
	
	for(i = 0; i < cnt; i++)
	{	
		ipt = &pointsData[i];
		//ipt->x <<= ipt->pyrLevel;
		//ipt->y <<= ipt->pyrLevel;
		
		const unsigned char *image = pyr_img[ipt->pyrLevel];
		
		fastRotationEstimation(image, ipt);
	}
	
	return cnt;
}


int supressNeighborPoints(IPoint *pointsData, const int cnt, const double maxDist)
{
	int i, j;		
	int neighb_map[cnt];
	
	qsort(pointsData, cnt, sizeof(IPoint), compareIPoint);
	
	for(i = 0; i < cnt; i++) neighb_map[i] = -1;
	
	for(i = 0; i < cnt; i++)
	{
		if(neighb_map[i] > -1) continue;
		
		const IPoint *p1 = &*(pointsData+i);
		for(j = i+1; j < cnt; j++)
		{
			if(neighb_map[j] > -1) continue;
			
			const IPoint *p2 = &*(pointsData+j);
			const int dx = p1->x - p2->x;
			const int dy = p1->y - p2->y;
			if( (dx*dx+dy*dy) < maxDist )
			{
				neighb_map[j] = i;
			}
		}
	}
	
	for(i = 0, j = 0; i < cnt; i++)
	{
		if(neighb_map[i] == -1)
		{
			memcpy( pointsData+j, pointsData+i, sizeof(IPoint) );
			j++;
		}
	}
	
	return j;
}

void fastRotationEstimation(const unsigned char *image, IPoint *ip)
{
	int i;
	double dx = 0.0;
	double dy = 0.0;
	
	const int pos = ip->pos;
	const double centrepx = (double)image[pos];
	
	const int stride = pyr_img1_width >> ip->pyrLevel;
	const double *ring_x = fast_ring_x;
	const double *ring_y = fast_ring_y;

	for (i = 0; i < 16; i++)
	{
		const double diff = (double)image[ (indY[i] * stride + indX[i]) + pos ] - centrepx;
		dx += diff * (*(ring_x++));
		dy += diff * (*(ring_y++));
	}

	ip->orientation = getCoterminalAngle( fast_atan2(dy, dx) );//ANGLE(dx, dy);
}

double FindShiTomasiScoreAtPoint(const unsigned char *image, const int stride, const xy *irCenter)
{
	double dXX = 0.0;
	double dYY = 0.0;
	double dXY = 0.0;
	
	const int cpos = irCenter->y * stride + irCenter->x;
	const int *win = shiTomasiWin;
	int len = 196 / 4;
	
	while(-- len > -1)
	{
		const double __dx0 = (double)image[ *(win++) + cpos ];
		const double __dx1 = (double)image[ *(win++) + cpos ];
		const double __dy0 = (double)image[ *(win++) + cpos ];
		const double __dy1 = (double)image[ *(win++) + cpos ];
		
		
		const double dx = __dx0 - __dx1;
		const double dy = __dy0 - __dy1;
		
		dXX += dx*dx;
		dYY += dy*dy;
		dXY += dx*dy;
	}
	
	const double nPixels = 0.002551020408163265;//1.0 / (2.0 * (double)((nx+1) * (ny+1)));
	dXX = dXX * nPixels;
	dYY = dYY * nPixels;
	dXY = dXY * nPixels;

	// Find and return smaller eigenvalue:
	return (double)0.5 * (dXX + dYY - fast_sqrt( (dXX + dYY) * (dXX + dYY) - (double)4.0 * (dXX * dYY - dXY * dXY) ));
}

void calcShiTomasiWin(const int stride, const int nHalfBoxSize)
{
	int stx = 0 - nHalfBoxSize;
	int sty = 0 - nHalfBoxSize;
	int enx = 0 + nHalfBoxSize;
	int eny = 0 + nHalfBoxSize;
	int y, x, j = 0;
	
	for(y = sty; y <= eny; y++)
	{
		for(x = stx; x <= enx; x++)
		{
			
			shiTomasiWin[j++] = y*stride + x + 1;
			shiTomasiWin[j++] = y*stride + x - 1;
			shiTomasiWin[j++] = y*stride + x + stride;
			shiTomasiWin[j++] = y*stride + x - stride;
		}
	}
}