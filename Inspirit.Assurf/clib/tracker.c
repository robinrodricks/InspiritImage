
#include "lkpyramid.h"

void getPointSample(const unsigned char *image, const int stride, IPoint *ipt);
void getPointSampleRotated(const unsigned char *image, const int stride, IPoint *ipt);
int findPrevFrameMatches(IPoint *pointsData, const int cnt, IPointMatch *prevFrameMatches, const int prevCnt);
int getClosestPoint(const IPoint *pointsData, const int cnt, const double x, const double y, const double maxDist);
double compareNCC(IPoint *a, IPoint *b);
int bestNCCMatch(IPoint *k, IPoint *pointsData, int *indices, const int cnt);
int searchClosestPoints(const IPoint *pointsData, const int cnt, const int x, const int y, const int maxMotion, int *indices);


int findPrevFrameMatches(IPoint *pointsData, const int cnt, IPointMatch *prevFrameMatches, const int prevCnt)
{
	const int nMaxSSDPerPixel = 500;
	const int mnMaxSSD = 16 * 16 * nMaxSSDPerPixel;
	const int distThresh = 12 * 12;
	
	int i, j, k;
	
	int dist_map[prevCnt];
	int sim_map[prevCnt];
	int sampled;
	int dist_b;
	int best_ind;
	IPoint *pt_prev, *pt_new;
	
	for(i = 0; i < prevCnt; i++)
	{
		dist_map[i] = mnMaxSSD + 1;
		sim_map[i] = -1;
	}
	
	for(i = 0; i < cnt; i++)
	{
		pt_new = &pointsData[i];
		
		const unsigned char *curr_img = pyr_img_blur[pt_new->pyrLevel];
		
		sampled = 0;
		dist_b = mnMaxSSD + 1;
		best_ind = -1;
		
			
		for(j = 0; j < prevCnt; j++)
		{
			pt_prev = prevFrameMatches[j].first;
			
			const unsigned char *prev_img = pyr_img_prev[pt_prev->pyrLevel];
			
			if((iSquare(pt_prev->x+pt_prev->dx - pt_new->x) + iSquare(pt_prev->y+pt_prev->dy - pt_new->y)) < distThresh)
			{
				if(sampled == 0)
				{
					getPointSample(curr_img, pyr_img1_width>>pt_new->pyrLevel, pt_new);
					//getPointSampleRotated(curr_img, pyr_img1_width>>pt_new->pyrLevel, pt_new);
					pt_new->sampled = sampled = 1;
				}
				if(pt_prev->sampled == 0)
				{
					getPointSample(prev_img, pyr_img1_width>>pt_prev->pyrLevel, pt_prev);
					//getPointSampleRotated(prev_img, pyr_img1_width>>pt_prev->pyrLevel, pt_prev);
					pt_prev->sampled = 1;
				}
				
				int dist = 0.0;
				const unsigned char *desc1 = &*pt_prev->sample;
				const unsigned char *desc2 = &*pt_new->sample;
				k = 256;
				while( --k > -1 )
				{
					dist += (*(desc1++)) * (*(desc2++));
				}
				
				const int SA = pt_prev->sum;
				const int SB = pt_new->sum;

				dist = ((2.0*SA*SB - SA*SA - SB*SB)/256.0 + pt_new->sqsum + pt_prev->sqsum - 2.0*dist);

				if(dist_b > dist && dist_map[j] > dist)
				{
					best_ind = j;
					dist_map[j] = dist_b = dist;
				}
			}
		}
		
		if(best_ind > -1) sim_map[best_ind] = i;
	}
	
	k = 0;
	
	for(j = 0; j < prevCnt; j++)
	{
		i = sim_map[j];
		if(i > -1)
		{
			const IPoint *prev = prevFrameMatches[j].first;
			const int px = prev->x;
			const int py = prev->y;
			
			memcpy(prevFrameMatches + k, prevFrameMatches + j, sizeof(IPointMatch));
			
			pt_new = &pointsData[i];
			pt_new->dx = (pt_new->x - px);
			pt_new->dy = (pt_new->y - py);
			pt_new->matched = 1;
			
			prevFrameMatches[k].first = pt_new;
			//prevFrameMatches[j].first = pt_new;
			//prevFrameMatches[j].tracked++;
			
			k++;
		}
	}
	
	return k;
}

int ncc_calls;

int NCCLKTrack(IPoint *pointsData, int cnt, IPointMatch *prevFrameMatches, const int prevCnt)
{
	const int MAX_TRACK_FT = 256;
	const int maxMotion = 10 * 10;
	int i, j;
	int ncnt = cnt;
	int nft = 0;
	
	ncc_calls = 0;
	int nmatches = 0;
	int nsaved = 0;
	
	IPointMatch *match;
	
	// NCC Frame to Frame
	int indices[cnt];
	for (i = 0; i < prevCnt; i++) 
	{
		match = &prevFrameMatches[i];
		IPoint *k = match->first;
		
		int kn = searchClosestPoints(pointsData, cnt, k->x + k->dx, k->y + k->dy, maxMotion, &*indices);
		int match_idx = bestNCCMatch(k, pointsData, indices, kn);
		if(match_idx > -1) 
		{
			IPoint *mpt = &pointsData[match_idx];
			if (mpt->matched == 0)
			{
				mpt->dx = (mpt->x - k->x);
				mpt->dy = (mpt->y - k->y);
				mpt->matched = 1;
				
				match->prev = k;
				match->first = mpt;
				match->tracked++;
				nmatches++;
			} else {
				for(j = 0; j < i; j++)
				{
					if(prevFrameMatches[j].first->matched && prevFrameMatches[j].first->index == mpt->index)
					{
						prevFrameMatches[j].first = prevFrameMatches[j].prev;
						break;
					}
				}
				nmatches--;
				mpt->matched = 0;
			}
		}
	}
	
	// LK Optical Flow
	CvPoint2D32f prev_ft[MAX_TRACK_FT];
	int prev_kpt[MAX_TRACK_FT];
	char lk_status[MAX_TRACK_FT];
	CvPoint2D32f curr_ft[MAX_TRACK_FT];
	
	// try to track lost keypoints using template matching
	for (i = 0; i < prevCnt; i++) 
	{
		match = &prevFrameMatches[i];
		IPoint *k = match->first;
		
		if ( k->matched == 0 && match->tracked > 2 && match->lk_tracked <= match->tracked>>1 ) 
		{
			prev_ft[nft].x = k->x;
			prev_ft[nft].y = k->y;
			// prediction
			curr_ft[nft].x = (k->x + k->dx);
			curr_ft[nft].y = (k->y + k->dy); 
			prev_kpt[nft] = i;
			nft++;
			
			if (nft>=MAX_TRACK_FT) break;
		} 
	}
	
	//unsigned char *my_pyr_img[3] = { pyr_img1, pyr_img_blur2, pyr_img_blur3 };
	const int sizeW[3] = {pyr_img1_width, pyr_img2_width, pyr_img3_width};
	const int sizeH[3] = {pyr_img1_height, pyr_img2_height, pyr_img3_height};
	
	myCalcOpticalFlowPyrLK(pyr_img_prev, /*pyr_img*/pyr_img_blur, sizeW, sizeH, 
							&*prev_ft, &*curr_ft, nft, 16, &*lk_status,
							0, CV_LKFLOW_INITIAL_GUESSES);
	
	for (i = 0; i < nft; i++) 
	{
		if (lk_status[i] == 0) continue;
		
		match = &prevFrameMatches[ prev_kpt[i] ];
		IPoint *prev = match->first;
		const int px = prev->x;
		const int py = prev->y;

		int closest_idx = getClosestPoint(pointsData, cnt, curr_ft[i].x, curr_ft[i].y, 4.0);
		if (closest_idx > -1) 
		{
			IPoint *pt = &pointsData[closest_idx];
			if (pt->matched == 0) 
			{
				//set_match(prev_kpt[i], pt);
				
				pt->dx = (pt->x - px);
				pt->dy = (pt->y - py);
				pt->matched = 1;
				
				match->prev = prev;
				match->first = pt;
				match->lk_tracked++;
				match->tracked++;
				nsaved++;
			}
		}
		else
		{
			// NCC TRY
			IPoint *newkpt = &pointsData[ncnt];
			
			newkpt->x = curr_ft[i].x;
			newkpt->y = curr_ft[i].y;
			newkpt->dx = newkpt->dy = 0;
			newkpt->pos = (newkpt->y>>prev->pyrLevel) * sizeW[prev->pyrLevel] + (newkpt->x>>prev->pyrLevel);
			newkpt->pyrLevel = prev->pyrLevel;
			newkpt->score = prev->score;
			newkpt->orientation = prev->orientation;
			newkpt->sampled = 0;
			newkpt->matched = 0;
			if (compareNCC(prev, newkpt) > 0.88)
			{
				//set_match(prev_kpt[i], newkpt);
				
				newkpt->dx = (newkpt->x - px);
				newkpt->dy = (newkpt->y - py);
				newkpt->matched = 1;
				
				match->prev = prev;
				match->first = newkpt;
				match->lk_tracked++;
				match->tracked++;
				ncnt++;
				nsaved++;
			}
		}
	}
	
	int mcnt = 0;
	for (i = 0; i < prevCnt; i++) 
	{
		IPoint *k = prevFrameMatches[i].first;
		
		if (k->matched == 1) 
		{
			memcpy(prevFrameMatches + mcnt, prevFrameMatches + i, sizeof(IPointMatch));
			mcnt++;
		} 
	}
	
	/*
	cout << ncc_calls << " calls to cmp_ncc, for " << cnt << " points in frame. Avg: " << ncc_calls / cnt;
	
	cout << nmatches << " ncc matches";
	cout << ", tracked " << nft << " features with LK. Saved " << nsaved << " tracks. ";
	cout << "tot: " << nmatches + nsaved << " features followed. " << (100.0*nsaved/(nsaved+nmatches)) << "% LK tracked.\n";
	*/
	result_b[0] = ncc_calls;
	result_b[1] = nmatches;
	result_b[2] = nft;
	result_b[3] = nsaved;
	
	return mcnt;
}

int getClosestPoint(const IPoint *pointsData, const int cnt, const double x, const double y, const double maxDist)
{
	int i;
	double best_dist = maxDist;
	int best_idx = -1;
	
	for(i = 0; i < cnt; i++)
	{		
		const IPoint *p1 = &*(pointsData+i);
		const double dx = (double)p1->x - x;
		const double dy = (double)p1->y - y;
		const double dist = (dx*dx+dy*dy);
		if( dist < best_dist )
		{
			best_dist = dist;
			best_idx = i;
		}
	}
	
	return best_idx;
}

int searchClosestPoints(const IPoint *pointsData, const int cnt, const int x, const int y, const int maxMotion, int *indices)
{
	int i, n = 0;
	
	for(i = 0; i < cnt; i++)
	{		
		const IPoint *p1 = &*(pointsData+i);
		const int dx = p1->x - x;
		const int dy = p1->y - y;
		const int dist = (dx*dx+dy*dy);
		if( dist < maxMotion )
		{
			indices[n++] = i;
		}
	}
	
	return n;
}

double compareNCC(IPoint *a, IPoint *b)
{
	ncc_calls++;
	
	if(b->sampled == 0)
	{
		getPointSample(/*pyr_img*/pyr_img_blur[b->pyrLevel], pyr_img1_width>>b->pyrLevel, b);
		//getPointSampleRotated(/*pyr_img*/pyr_img_blur[b->pyrLevel], pyr_img1_width>>b->pyrLevel, b);
		b->sampled = 1;
	}
	if(a->sampled == 0)
	{
		getPointSample(pyr_img_prev[a->pyrLevel], pyr_img1_width>>a->pyrLevel, a);
		//getPointSampleRotated(pyr_img_prev[a->pyrLevel], pyr_img1_width>>a->pyrLevel, a);
		a->sampled = 1;
	}
	
	int ma = a->mean;
	int mb = b->mean;
	int sum = 0;
	double norm = a->stdev * b->stdev;
	
	const unsigned char *desc1 = &*a->sample;
	const unsigned char *desc2 = &*b->sample;
	int k = 256;
	while( --k > -1 )
	{
		sum += (*(desc1++)-ma) * (*(desc2++)-mb);
	}

	return (double)sum / norm;
}

int bestNCCMatch(IPoint *k, IPoint *pointsData, int *indices, const int cnt)
{
	int i;
	double best_corr = 0.88;
	int best_idx = -1;

	for (i = 0; i < cnt; i++)
	{
		IPoint *pt = &pointsData[ indices[i] ];
		double ncc = compareNCC(k, pt);
		if (ncc > best_corr)
		{
			best_corr = ncc;
			best_idx = indices[i];
			if (best_corr > 0.9) break;
		}
	}

	return best_idx;
}

void getPointSample(const unsigned char *image, const int stride, IPoint *ipt)
{
	int k, m;
	const int SAMPLE_WIN = 16;
	const int SAMPLE_SIZE = SAMPLE_WIN * SAMPLE_WIN;
	const int halfWin = SAMPLE_WIN / 2;
	
	int ss = 0;
	int sq = 0;
	
	unsigned char *sample = &*ipt->sample;
	const int cpos = ipt->pos;
	
	for(k = 0; k < SAMPLE_WIN; k++)
	{
		for(m = 0; m < SAMPLE_WIN; m++)
		{
			const int ox = m - halfWin;
			const int oy = k - halfWin;
			
			const unsigned char val = image[ cpos + ( oy * stride + ox ) ];
			ss += val;
			sq += (int)(val * val);
				
			*(sample++) = val;
		}
	}
	
	ipt->sum = ss;
	ipt->sqsum = sq;
	ipt->mean = ss / SAMPLE_SIZE;
	ipt->stdev = fast_sqrt( sq - (SAMPLE_SIZE * ipt->mean * ipt->mean) );
}

void getPointSampleRotated(const unsigned char *image, const int stride, IPoint *ipt)
{	
	const int SAMPLE_WIN = 16;
	const int SAMPLE_SIZE = SAMPLE_WIN * SAMPLE_WIN;
	const int halfWin = SAMPLE_WIN / 2;
	int m, k;
	
	double si, co;
	sin_cos(ipt->orientation, &si, &co);
	
	const double origx = ipt->x>>ipt->pyrLevel;
	const double origy = ipt->y>>ipt->pyrLevel;
	
	int ss = 0;
	int sq = 0;
	
	unsigned char *sample = &*ipt->sample;
	
	for(k = 0; k < SAMPLE_WIN; k++)
	{
		for(m = 0; m < SAMPLE_WIN; m++)
		{
			const int ox = m - halfWin;
			const int oy = k - halfWin;
			
			const double rotX = (co * (double)ox - si * (double)oy) + origx;
			const double rotY = (si * (double)ox + co * (double)oy) + origy;
			const unsigned char val = image[ dRound(rotX) + dRound(rotY) * stride ];
			//const unsigned char val = (unsigned char)bilinear_interpolation(image, stride, rotX, rotY);
			
			ss += val;
			sq += (int)(val * val);
			
			*(sample++) = val;
		}
	}
	ipt->sum = ss;
	ipt->sqsum = sq;
	ipt->mean = ss / SAMPLE_SIZE;
	ipt->stdev = fast_sqrt( sq - (SAMPLE_SIZE * ipt->mean * ipt->mean) );
}