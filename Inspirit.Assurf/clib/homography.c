

int findHomography(const int np, register double *obj, register double *img, register double *mat);
int matchesHomography(IPointMatch *matches, const int np, register double *mat);


#include "svdcmp.c"
#include "cmphomo.c"

void locatePlanarObject(IPointMatch *matches, int *matchesCount, double *homography)
{
	int np = *matchesCount;
	if(np > 4)
	{
		qsort(matches, np, sizeof(IPointMatch), compareIPointMatch);
		
		int i, j;
		double u_v_up_vp[np*4];
		double normalized_u_v_up_vp[np*4];
		double T1[9];
		double T2inv[9];
		double tmp[9];
		double bestH[9];

		for(i = 0, j = 0; i < np; i++)
		{
			IPointMatch *mtch = &matches[i];
			mtch->first->matched = 0;
			mtch->wasGood = 0;
			u_v_up_vp[j++] = (double)mtch->second->x;
			u_v_up_vp[j++] = (double)mtch->second->y;
			u_v_up_vp[j++] = (double)mtch->first->x;
			u_v_up_vp[j++] = (double)mtch->first->y;
		}
		
		normalizePoints(np, u_v_up_vp, &*normalized_u_v_up_vp, &*T1, &*T2inv);
		
		const double logProb = log( (double)1.0 - PROBABILITY_REQUIRED );
		
		int N = np > 10 ? 100 : 1000;
		int number_of_inliers = 0;
		int sample_count = 0;
		int prosac_correspondences = np >= 10 ? 10 : 5;
		int current_inliers[np];
		int inliers[np];

		int n1, n2, n3, n4;
		
		while (N > sample_count && number_of_inliers < 20)
		{
			get_4_prosac_indices(prosac_correspondences, &n1, &n2, &n3, &n4);

			// This incrementing strategy is naive and simple but works just fine most of the time.
			if(prosac_correspondences < np) {
				++prosac_correspondences;
			}

			homography_from_4corresp(homography,
							   normalized_u_v_up_vp[n1], normalized_u_v_up_vp[n1 + 1], normalized_u_v_up_vp[n1 + 2], normalized_u_v_up_vp[n1 + 3],
							   normalized_u_v_up_vp[n2], normalized_u_v_up_vp[n2 + 1], normalized_u_v_up_vp[n2 + 2], normalized_u_v_up_vp[n2 + 3],
							   normalized_u_v_up_vp[n3], normalized_u_v_up_vp[n3 + 1], normalized_u_v_up_vp[n3 + 2], normalized_u_v_up_vp[n3 + 3],
							   normalized_u_v_up_vp[n4], normalized_u_v_up_vp[n4 + 1], normalized_u_v_up_vp[n4 + 2], normalized_u_v_up_vp[n4 + 3]);

			//denormalize Homography
			MultiplyMat(T2inv, homography, &*tmp, 3, 3, 3, 3);
			MultiplyMat(tmp, T1, homography, 3, 3, 3, 3);			
			set_bottom_right_coefficient_to_one(homography);
			
			if (nice_homography(homography)) 
			{
				int current_number_of_inliers = compute_inliers(matches, np, homography, &*current_inliers, INLIER_THRESHOLD_SQ);

				if (current_number_of_inliers > number_of_inliers) 
				{					
					double inv_epsilon = (double)1.0 - ((double)1.0 - ((double)current_number_of_inliers) / ((double)np));
					int temp = (int)(logProb / log( (double)1.0 - (inv_epsilon * inv_epsilon * inv_epsilon * inv_epsilon) ) );
					if(temp > 0 && temp < N){
						N = temp;
					}

					number_of_inliers = current_number_of_inliers;
					memcpy(&*inliers, &*current_inliers, number_of_inliers * sizeof(int));
					memcpy(&*bestH, homography, 9 * sizeof(double));
				}
			}
			
			sample_count++;
		}
		
		/*
		double pts1[number_of_inliers*2];
		double pts2[number_of_inliers*2];
		
		for(i = 0; i < number_of_inliers; i++)
		{
			j = inliers[i] * 4;
			pts2[i*2] = normalized_u_v_up_vp[ j++ ];
			pts2[i*2+1] = normalized_u_v_up_vp[ j++ ];
			pts1[i*2] = normalized_u_v_up_vp[ j++ ];
			pts1[i*2+1] = normalized_u_v_up_vp[ j++ ];
		}
		
		findHomography(number_of_inliers, pts1, pts2, homography);
		MultiplyMat(T2inv, homography, &*tmp, 3, 3, 3, 3);
		MultiplyMat(tmp, T1, homography, 3, 3, 3, 3);		
		set_bottom_right_coefficient_to_one(homography);
		
		number_of_inliers = compute_inliers(matches, np, homography, &*inliers, INLIER_THRESHOLD_SQ);
		*/
		memcpy(homography, &*bestH, 9 * sizeof(double));
		
		for(i = 0; i < number_of_inliers; i++)
		{
			memcpy(matches + i, matches + inliers[i], sizeof(IPointMatch));
			matches[i].first->matched = 1;
			matches[i].wasGood = 1;
		}
		
		*matchesCount = number_of_inliers;
	}
}

const int homo_16_idx[16][4] = {
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

void locatePlanarObject16(IPointMatch *matches, int *matchesCount, double *homography)
{
	int np = *matchesCount;
	if(np > 4)
	{
		qsort(matches, np, sizeof(IPointMatch), compareIPointMatch);
		
		int i, j;
		double u_v_up_vp[np*4];
		double normalized_u_v_up_vp[np*4];
		double T1[9];
		double T2inv[9];
		double tmp[9];
		double bestH[9];
		
		int N;
		int number_of_inliers = 0;
		int sample_count = 0;
		int current_inliers[np];
		int inliers[np];
		
		int n1, n2, n3, n4;
		
		for(i = 0, j = 0; i < np; i++)
		{
			IPointMatch *mtch = &matches[i];
			mtch->first->matched = 0;
			u_v_up_vp[j++] = (double)mtch->second->x;
			u_v_up_vp[j++] = (double)mtch->second->y;
			u_v_up_vp[j++] = (double)mtch->first->x;
			u_v_up_vp[j++] = (double)mtch->first->y;
		}
		
		normalizePoints(np, u_v_up_vp, &*normalized_u_v_up_vp, &*T1, &*T2inv);
		
		if(np >= 16)
		{
		
			int result_pt[8];
			perpendicularRegressionIdx(matches, np, &*result_pt);
			
			N = 16;
			
			while (N > sample_count)
			{
				n1 = result_pt[ homo_16_idx[sample_count][0] ]<<2;
				n2 = result_pt[ homo_16_idx[sample_count][1] ]<<2;
				n3 = result_pt[ homo_16_idx[sample_count][2] ]<<2;
				n4 = result_pt[ homo_16_idx[sample_count][3] ]<<2;

				homography_from_4corresp(homography,
							   normalized_u_v_up_vp[n1], normalized_u_v_up_vp[n1 + 1], normalized_u_v_up_vp[n1 + 2], normalized_u_v_up_vp[n1 + 3],
							   normalized_u_v_up_vp[n2], normalized_u_v_up_vp[n2 + 1], normalized_u_v_up_vp[n2 + 2], normalized_u_v_up_vp[n2 + 3],
							   normalized_u_v_up_vp[n3], normalized_u_v_up_vp[n3 + 1], normalized_u_v_up_vp[n3 + 2], normalized_u_v_up_vp[n3 + 3],
							   normalized_u_v_up_vp[n4], normalized_u_v_up_vp[n4 + 1], normalized_u_v_up_vp[n4 + 2], normalized_u_v_up_vp[n4 + 3]);

				MultiplyMat(T2inv, homography, &*tmp, 3, 3, 3, 3);
				MultiplyMat(tmp, T1, homography, 3, 3, 3, 3);			
				set_bottom_right_coefficient_to_one(homography);
				
				if (nice_homography(homography)) 
				{
					int current_number_of_inliers = compute_inliers(matches, np, homography, &*current_inliers, INLIER_THRESHOLD_SQ);

					if (current_number_of_inliers > number_of_inliers) 
					{
						number_of_inliers = current_number_of_inliers;
						memcpy(&*inliers, &*current_inliers, number_of_inliers * sizeof(int));
						memcpy(&*bestH, homography, 9 * sizeof(double));
					}
				}
				
				sample_count++;
			}
		}
		else
		{			
			const double logProb = log( (double)1.0 - PROBABILITY_REQUIRED );
			N = 1000;
			int prosac_correspondences = 5;
			
			while (N > sample_count)
			{
				get_4_prosac_indices(prosac_correspondences, &n1, &n2, &n3, &n4);

				// This incrementing strategy is naive and simple but works just fine most of the time.
				if(prosac_correspondences < np) {
					++prosac_correspondences;
				}

				homography_from_4corresp(homography,
							   normalized_u_v_up_vp[n1], normalized_u_v_up_vp[n1 + 1], normalized_u_v_up_vp[n1 + 2], normalized_u_v_up_vp[n1 + 3],
							   normalized_u_v_up_vp[n2], normalized_u_v_up_vp[n2 + 1], normalized_u_v_up_vp[n2 + 2], normalized_u_v_up_vp[n2 + 3],
							   normalized_u_v_up_vp[n3], normalized_u_v_up_vp[n3 + 1], normalized_u_v_up_vp[n3 + 2], normalized_u_v_up_vp[n3 + 3],
							   normalized_u_v_up_vp[n4], normalized_u_v_up_vp[n4 + 1], normalized_u_v_up_vp[n4 + 2], normalized_u_v_up_vp[n4 + 3]);

				MultiplyMat(T2inv, homography, &*tmp, 3, 3, 3, 3);
				MultiplyMat(tmp, T1, homography, 3, 3, 3, 3);			
				set_bottom_right_coefficient_to_one(homography);
				
				if (nice_homography(homography)) 
				{
					int current_number_of_inliers = compute_inliers(matches, np, homography, &*current_inliers, INLIER_THRESHOLD_SQ);

					if (current_number_of_inliers > number_of_inliers) 
					{					
						double inv_epsilon = (double)1.0 - ((double)1.0 - ((double)current_number_of_inliers) / ((double)np));
						int temp = (int)(logProb / log( (double)1.0 - (inv_epsilon * inv_epsilon * inv_epsilon * inv_epsilon) ) );
						if(temp > 0 && temp < N){
							N = temp;
						}

						number_of_inliers = current_number_of_inliers;
						memcpy(&*inliers, &*current_inliers, number_of_inliers * sizeof(int));
						memcpy(&*bestH, homography, 9 * sizeof(double));
					}
				}
				
				sample_count++;
			}
		}
		/*
		// compute homography using best inlier set
		double pts1[number_of_inliers*2];
		double pts2[number_of_inliers*2];
		
		for(i = 0; i < number_of_inliers; i++)
		{
			j = inliers[i] * 4;
			pts2[i*2] = normalized_u_v_up_vp[ j++ ];
			pts2[i*2+1] = normalized_u_v_up_vp[ j++ ];
			pts1[i*2] = normalized_u_v_up_vp[ j++ ];
			pts1[i*2+1] = normalized_u_v_up_vp[ j++ ];
		}
		
		findHomography(number_of_inliers, pts1, pts2, homography);
		
		MultiplyMat(T2inv, homography, &*tmp, 3, 3, 3, 3);
		MultiplyMat(tmp, T1, homography, 3, 3, 3, 3);		
		set_bottom_right_coefficient_to_one(homography);
		
		number_of_inliers = compute_inliers(matches, np, homography, &*inliers, INLIER_THRESHOLD_SQ);
		*/
		memcpy(homography, &*bestH, 9 * sizeof(double));
		for(i = 0; i < number_of_inliers; i++)
		{
			memcpy(matches + i, matches + inliers[i], sizeof(IPointMatch));
			matches[i].first->matched = 1;
		}
		
		*matchesCount = number_of_inliers;
	}
}

void directHomography(IPointMatch *matches, int *matchesCount, double *homography)
{
	int np = *matchesCount;
	
	int i, j;
	double u_v_up_vp[np*4];
	double normalized_u_v_up_vp[np*4];
	double T1[9];
	double T2inv[9];
	double tmp[9];
	
	int number_of_inliers = np;
	int inliers[np];
	
	for(i = 0, j = 0; i < np; i++)
	{
		const IPointMatch *mtch = &matches[i];
		u_v_up_vp[j++] = (double)mtch->second->x;
		u_v_up_vp[j++] = (double)mtch->second->y;
		u_v_up_vp[j++] = (double)mtch->first->x;
		u_v_up_vp[j++] = (double)mtch->first->y;
	}
	
	normalizePoints(np, u_v_up_vp, &*normalized_u_v_up_vp, &*T1, &*T2inv);
	
	double pts1[number_of_inliers*2];
	double pts2[number_of_inliers*2];
	
	for(i = 0, j = 0; i < number_of_inliers; i++)
	{
		pts2[i*2] = normalized_u_v_up_vp[ j++ ];
		pts2[i*2+1] = normalized_u_v_up_vp[ j++ ];
		pts1[i*2] = normalized_u_v_up_vp[ j++ ];
		pts1[i*2+1] = normalized_u_v_up_vp[ j++ ];
	}
	
	findHomography(number_of_inliers, pts1, pts2, homography);
	
	MultiplyMat(T2inv, homography, &*tmp, 3, 3, 3, 3);
	MultiplyMat(tmp, T1, homography, 3, 3, 3, 3);		
	set_bottom_right_coefficient_to_one(homography);
	
	number_of_inliers = compute_inliers(matches, np, homography, &*inliers, INLIER_THRESHOLD_SQ);
	
	for(i = 0; i < number_of_inliers; i++)
	{
		memcpy(matches + i, matches + inliers[i], sizeof(IPointMatch));
	}
	
	*matchesCount = number_of_inliers;
}

int findHomography(const int np, register double *pts1, register double *pts2, register double *mat)
{
	const int np2 = np << 1;

	double a[np2*8], b[np2], temp[np2*8];

	int i, j, ii, jj;
	double sx, sy, dx, dy;

	for( i = 0, j = np; i < np; i++, j++ )
	{
		dx = *(pts1++);
		dy = *(pts1++);
		sx = *(pts2++);
		sy = *(pts2++);
		
		ii = i<<3;
		jj = j<<3;
		
		a[ii] = a[jj+3] = sx;
		a[ii+1] = a[jj+4] = sy;
		a[ii+2] = a[jj+5] = 1;
		a[ii+3] = a[ii+4] = a[ii+5] =
		a[jj] = a[jj+1] = a[jj+2] = 0.0;
		a[ii+6] = -dx*sx;
		a[ii+7] = -dx*sy;
		a[jj+6] = -dy*sx;
		a[jj+7] = -dy*sy;
		b[i]    = dx;
		b[j]    = dy;
	}

	if(PseudoInverse(&*temp, a, np2, 8)){
		//return 1;
	}
	MultiplyMat(temp, b, &*mat, 8, np2, np2, 1);
	//mat[8] = 1;
	
	return 0;
}

int matchesHomography(IPointMatch *matches, const int np, register double *mat)
{
	const int np2 = np << 1;

	double a[np2*8], b[np2], temp[np2*8];

	int i, j, ii, jj;
	double sx, sy, dx, dy;

	for( i = 0, j = np; i < np; i++, j++ )
	{
		dx = matches[i].first->x;
		dy = matches[i].first->y;
		sx = matches[i].second->x;
		sy = matches[i].second->y;
		
		ii = i<<3;
		jj = j<<3;
		
		a[ii] = a[jj+3] = sx;
		a[ii+1] = a[jj+4] = sy;
		a[ii+2] = a[jj+5] = 1;
		a[ii+3] = a[ii+4] = a[ii+5] =
		a[jj] = a[jj+1] = a[jj+2] = 0.0;
		a[ii+6] = -dx*sx;
		a[ii+7] = -dx*sy;
		a[jj+6] = -dy*sx;
		a[jj+7] = -dy*sy;
		b[i]    = dx;
		b[j]    = dy;
	}

	if(PseudoInverse(&*temp, a, np2, 8)){
		return 1;
	}
	MultiplyMat(temp, b, &*mat, 8, np2, np2, 1);
	mat[8] = 1;
	
	return 0;
}