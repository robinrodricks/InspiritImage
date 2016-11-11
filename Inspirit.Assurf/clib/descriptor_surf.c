/*
void getPointDescriptor_SURF(const int x, const int y, const double scale, const double co, const double si, 
							register const double *integralData, double *start, double *end)
{
	int ind11, ind12, cc1, cc2;
	int sample_x, sample_y, k, l;
	int i, ix=0, j, jx=0, xs=0, ys=0;
	double scale25, dx, dy, mdx, mdy;
	double gauss_s1=0, gauss_s2=0;
	double rx=0.0, ry=0.0, rrx=0.0, rry=0.0, len=0.0;
	double cx = -0.5, cy = 0.0; //Subregion centers for the 4x4 IGAUSSIAN weighting

	scale25 = 2.5 * scale;
	int s2 = (dRound(scale)<<1);
	int s22 = (s2>>1);

	register double *data_ptr = start;

	//Calculate descriptor for this interest point
	for( i = -8; i < 12; i+=9 )
	{
		i = i-4;

		cx += 1.0;
		cy = -0.5;

		for( j = -8; j < 12; j+=9 )
		{
			dx=dy=mdx=mdy=0;
			cy += 1.0;

			j = j - 4;

			ix = i + 5;
			jx = j + 5;

			xs = dRound(x + ( -jx*scale*si + ix*scale*co ));
			ys = dRound(y + ( jx*scale*co + ix*scale*si ));

			for (k = i; k < i + 9; ++k)
			{
				for (l = j; l < j + 9; ++l)
				{
					//Get coords of sample point on the rotated axis
					sample_x = dRound(x + (-l*scale*si + k*scale*co));
					sample_y = dRound(y + ( l*scale*co + k*scale*si));

					//Get the GAUSSIAN weighted x and y responses
					gauss_s1 = IGAUSSIAN(xs-sample_x, ys-sample_y, scale25);

					sample_x += iborder - 1;
					sample_y += iborder - 1;

					rx = FMAX(0.0, (*(integralData+(ind11=(sample_y - s22)*iwidth)+sample_x) - *(integralData+ind11+(cc2=sample_x+s22)) - *(integralData+(ind12=(sample_y - s22+s2)*iwidth)+sample_x) + *(integralData+ind12+cc2)))
					- 1*FMAX(0.0, (*(integralData+ind11+(cc1=sample_x-s22)) - *(integralData+ind11+sample_x) - *(integralData+ind12+cc1) + *(integralData+ind12+sample_x)));

					ry = FMAX(0.0, (*(integralData+(ind11=(sample_y)*iwidth)+(cc1=sample_x-s22)) - *(integralData+ind11+(cc2=cc1+s2)) - *(integralData+(ind12=(sample_y + s22)*iwidth)+cc1) + *(integralData+ind12+cc2)))
					- 1*FMAX(0.0, (*(integralData+(ind11=(sample_y - s22)*iwidth)+(cc1=sample_x-s22)) - *(integralData+ind11+(cc2=cc1+s2)) - *(integralData+(ind12=(sample_y)*iwidth)+cc1) + *(integralData+ind12+cc2)));


					//Get the IGAUSSIAN weighted x and y responses on rotated axis
					rrx = gauss_s1 * (-rx*si + ry*co);
					rry = gauss_s1 * (rx*co + ry*si);

					dx += rrx;
					dy += rry;
					mdx += fabs(rrx);
					mdy += fabs(rry);
				}
			}

			gauss_s2 = FGAUSSIAN(cx-2.0, cy-2.0, 1.5);

			*(data_ptr++) = dx*gauss_s2;
			*(data_ptr++) = dy*gauss_s2;
			*(data_ptr++) = mdx*gauss_s2;
			*(data_ptr++) = mdy*gauss_s2;

			len += (dx*dx + dy*dy + mdx*mdx + mdy*mdy) * gauss_s2*gauss_s2;

		}
	}

	//Convert to Unit Vector
	len = 1.0 / fast_sqrt(len);
	for(data_ptr = start; data_ptr < end;)
	{
		*(data_ptr++) *= len;
	}
}
*/

/* Gaussian used to weight descriptor samples */
static double gaussWeights20[20][20];

void calcGaussWeights(const int PATCH_SZ)
{
	const double DESC_SIGMA = 3.3;
	const double halfPath = (double)(PATCH_SZ-1)*0.5;
	int i, j;
    double c2 = 1./(DESC_SIGMA*DESC_SIGMA*2.0);
    double gs = 0.;
	for( i = 0; i < PATCH_SZ; i++ )
    {
        for( j = 0; j < PATCH_SZ; j++ )
        {
            double __dx = (double)j - halfPath;
			double __dy = (double)i - halfPath;
            double val = exp(-(__dx*__dx+__dy*__dy)*c2);
            gaussWeights20[i][j] = val;
            gs += val;
        }
    }
	gs = (double)1. / gs;
	for( i = 0; i < PATCH_SZ; i++ )
    {
        for( j = 0; j < PATCH_SZ; j++ )
        {
			gaussWeights20[i][j] *= gs;
		}
	}
}
void SURFDescriptors( IPoint *pointsData, const int N, double *descrData, const int extended )
{
	const int PATCH_SZ = 20;
	const int win_size = (int)((PATCH_SZ+1)*2.0);
	const int d = 5;
	const int descrSize = (extended == 1) ? 128 : 64;
	
	int i, j, k;

	IPoint *ipt;
	int x, y;
	double PATCH[PATCH_SZ+1][PATCH_SZ+1];
	double DX[PATCH_SZ][PATCH_SZ];
	double DY[PATCH_SZ][PATCH_SZ];
	unsigned char WIN[win_size*win_size];
	double __dx, __dy, _mdx, _mdy, len, mean;
	double __dx2, __dy2, _mdx2, _mdy2;
	double pixel_x, pixel_y, start_x, start_y;
	double sin_dir, cos_dir;
	
	register double *descr = descrData;
	double *data_end;
	
	const double win_offset = -((double)(win_size-1) * (double)0.5);
	//const double win_offset = -((double)(PATCH_SZ-1) * (double)0.5);
	//const double inv_len = (double)1.0 / (double)descrSize;
	
    for(k = 0; k < N; k++ )
    {		
		ipt = &pointsData[k];
		
		if(ipt->matched) continue;
		
		const unsigned char *image = pyr_img_blur[ipt->pyrLevel];
		const int width = pyr_img1_width>>ipt->pyrLevel;
		
		ipt->descriptor = descr;
		data_end = descr + descrSize;
		
		const double orientation = ipt->orientation;
		sin_cos(-orientation, &sin_dir, &cos_dir);

        start_x = (double)(ipt->x>>ipt->pyrLevel) + win_offset*cos_dir + win_offset*sin_dir;
        start_y = (double)(ipt->y>>ipt->pyrLevel) - win_offset*sin_dir + win_offset*cos_dir;
        for( i=0; i<win_size; i++, start_x+=sin_dir, start_y+=cos_dir )
        //for( i=0; i<PATCH_SZ+1; i++, start_x+=sin_dir, start_y+=cos_dir )
        {
            pixel_x = start_x;
            pixel_y = start_y;
			for( j=0; j<win_size; j++, pixel_x+=cos_dir, pixel_y-=sin_dir )
            //for( j=0; j<PATCH_SZ+1; j++, pixel_x+=cos_dir, pixel_y-=sin_dir )
            {
				WIN[i*win_size + j] = image[ dRound(pixel_x) + dRound(pixel_y) * width ];
				//PATCH[i][j] = bilinear_interpolation(image, width, pixel_x, pixel_y);
            }
        }
		
		// resize to patch
		int xd, yd;

		for(i = 0 ; i < win_size ; i += 2)
		{
			yd = i/2;
			const unsigned char *row = &*(WIN+(i * win_size));
			const unsigned char *nrow = &*(WIN+((i+1) * win_size));
			for(j = 0 ; j < win_size ; j += 2)
			{
				xd = j/2;
				PATCH[yd][xd] = (double)( row[j] + row[j+1] + nrow[j] + nrow[j+1] ) * (double)0.25;
			}        
		}
		
        // Calculate gradients in x and y with wavelets of size 2s
        for( i = 0; i < PATCH_SZ; i++ )
		{
            for( j = 0; j < PATCH_SZ; j++ )
            {
				const double dw = gaussWeights20[i][j];
				DX[i][j] = (double)(PATCH[i][j+1] - PATCH[i][j] + PATCH[i+1][j+1] - PATCH[i+1][j]) * dw;
				DY[i][j] = (double)(PATCH[i+1][j] - PATCH[i][j] + PATCH[i+1][j+1] - PATCH[i][j+1]) * dw;
            }
		}

        // Construct the descriptor
		
		len = 0.0;
		mean = 0.0;
		
		if(extended)
		{
			// 128 bin descriptor
			for(i = 0; i < 4; i++)
			{
				for(j = 0; j < 4; j++)
				{
					__dx = __dy = _mdx = _mdy = 0.0;
					__dx2 = __dy2 = _mdx2 = _mdy2 = 0.0;
					for(y = i*d; y < i*d+d; y++)
					{
						for(x = j*d; x < j*d+d; x++)
						{							
							const double tx = DX[y][x], ty = DY[y][x];
                            if( ty >= 0 )
                            {
                                __dx += tx;
                                _mdx += fabs(tx);
                            } else {
                                __dx2 += tx;
                                _mdx2 += fabs(tx);
                            }
                            if ( tx >= 0 )
                            {
                                __dy += ty;
                                _mdy += fabs(ty);
                            } else {
                                __dy2 += ty;
                                _mdy2 += fabs(ty);
                            }
						}
					}
					*(descr++) = __dx;
					*(descr++) = _mdx;
					*(descr++) = __dx2;
					*(descr++) = _mdx2;
					*(descr++) = __dy;
					*(descr++) = _mdy;
					*(descr++) = __dy2;
					*(descr++) = _mdy2;
					
					//mean += __dx + __dy + _mdx + _mdy + __dx2 + __dy2 + _mdx2 + _mdy2;
					len += (__dx*__dx + __dy*__dy + _mdx*_mdx + _mdy*_mdy + __dx2*__dx2 + __dy2*__dy2 + _mdx2*_mdx2 + _mdy2*_mdy2);
				}
			}
		}
		else
		{
			// 64 bin descriptor
			for(i = 0; i < 4; i++)
			{
				for(j = 0; j < 4; j++)
				{
					__dx = __dy = _mdx = _mdy = 0.0;
					for(y = i*d; y < i*d+d; y++)
					{
						for(x = j*d; x < j*d+d; x++)
						{
							__dx += DX[y][x];
							__dy += DY[y][x];
							_mdx += fabs(DX[y][x]);
							_mdy += fabs(DY[y][x]);
						}
					}
					*(descr++) = __dx;
					*(descr++) = __dy;
					*(descr++) = _mdx;
					*(descr++) = _mdy;
					
					//mean += __dx + __dy + _mdx + _mdy;
					len += (__dx*__dx + __dy*__dy + _mdx*_mdx + _mdy*_mdy);
				}
			}
		}
		
		// not that good as with 36 bin descriptor
		/*
		mean *= inv_len;
		len = (double)1.0 / ( fast_sqrt( len * inv_len - (mean * mean) ) + (double)1.0E-12 );
		for(descr -= descrSize; descr < data_end; descr++)
		{
			*(descr) = (*(descr) - mean) * len;
		}
		*/
		// original idea
		//Convert to Unit Vector
		len = (double)1.0 / ( fast_sqrt(len) + (double)1.0E-12 );
		for(descr -= descrSize; descr < data_end;)
		{
			*(descr++) *= len;
		}
    }
}

// modified FAST SURF descriptor
void fastSURFDescriptors( IPoint *pointsData, const int N, double *descrData )
{
	const int PATCH_SZ = 15;
	const int d = 5;
	
	int i, j, k;

	IPoint *ipt;
	int x, y;
	double PATCH[PATCH_SZ+1][PATCH_SZ+1];
	double DX[PATCH_SZ][PATCH_SZ];
	double DY[PATCH_SZ][PATCH_SZ];
	double __dx, __dy, _mdx, _mdy, len, mean;
	double pixel_x, pixel_y, start_x, start_y;
	double sin_dir, cos_dir;
	
	register double *descr = descrData;
	double *data_end;
	
	const double win_offset = -((double)(PATCH_SZ-1) * (double)0.5);
	const double inv_len = (double)1.0 / (double)36.0;
	
    for(k = 0; k < N; k++ )
    {		
		ipt = &pointsData[k];
		
		if(ipt->matched) continue;
		
		const unsigned char *image = pyr_img_blur[ipt->pyrLevel];
		const int width = pyr_img1_width>>ipt->pyrLevel;
		
		ipt->descriptor = descr;
		data_end = descr + 36;
		
		const double orientation = ipt->orientation;
		sin_cos(-orientation, &sin_dir, &cos_dir);

        start_x = (double)(ipt->x>>ipt->pyrLevel) + win_offset*cos_dir + win_offset*sin_dir;
        start_y = (double)(ipt->y>>ipt->pyrLevel) - win_offset*sin_dir + win_offset*cos_dir;
        for( i=0; i<PATCH_SZ+1; i++, start_x+=sin_dir, start_y+=cos_dir )
        {
            pixel_x = start_x;
            pixel_y = start_y;
            for( j=0; j<PATCH_SZ+1; j++, pixel_x+=cos_dir, pixel_y-=sin_dir )
            {
				PATCH[i][j] = bilinear_interpolation(image, width, pixel_x, pixel_y);
            }
        }

        // Calculate gradients in x and y with wavelets of size 2s
        for( i = 0; i < PATCH_SZ; i++ )
		{
            for( j = 0; j < PATCH_SZ; j++ )
            {
				//const double dw = gaussWeights20[i][j];
				DX[i][j] = (double)(PATCH[i][j+1] - PATCH[i][j] + PATCH[i+1][j+1] - PATCH[i+1][j]);// * dw;
				DY[i][j] = (double)(PATCH[i+1][j] - PATCH[i][j] + PATCH[i+1][j+1] - PATCH[i][j+1]);// * dw;
            }
		}

        // Construct the descriptor
		// 36-bin descriptor
		
		len = 0.0;
		mean = 0.0;
		
        for(i = 0; i < 3; i++)
		{
            for(j = 0; j < 3; j++)
            {
				__dx = __dy = _mdx = _mdy = 0.0;
                for(y = i*d; y < i*d+d; y++)
                {
                    for(x = j*d; x < j*d+d; x++)
                    {
						__dx += DX[y][x];
						__dy += DY[y][x];
						_mdx += fabs(DX[y][x]);
						_mdy += fabs(DY[y][x]);
                    }
                }
				*(descr++) = __dx;
				*(descr++) = __dy;
				*(descr++) = _mdx;
				*(descr++) = _mdy;
				
				mean += __dx + __dy + _mdx + _mdy;
				len += (__dx*__dx + __dy*__dy + _mdx*_mdx + _mdy*_mdy);
			}
		}
		
		// this works better cause of the smoother result
		
		mean *= inv_len;
		len = (double)1.0 / ( fast_sqrt( len * inv_len - (mean * mean) ) + (double)1.0E-12 );
		for(descr -= 36; descr < data_end; descr++)
		{
			*(descr) = (*(descr) - mean) * len;
		}
		
		/*
		//Convert to Unit Vector
		len = (double)1.0 / ( fast_sqrt(len) + (double)1.0E-12 );
		for(descr -= 36; descr < data_end;)
		{
			*(descr++) *= len;
		}*/
    }
}
/*
void calculateDescriptors_SURF64(IPoint *pointsData, const int count, const double *integral, double *descrData)
{
	int i;
	IPoint *ipt;
	
	for(i = 0; i < count; i++)
	{
		ipt = &pointsData[i];
		
		if(ipt->matched) continue;
		
		const double orientation = ipt->orientation;
		ipt->descriptor = descrData+(64 * i);
		
		double sn, cs;
		sin_cos(orientation, &sn, &cs);
		
		getPointDescriptor_SURF(ipt->x, ipt->y, 2.0, cs, sn, 
								integral, ipt->descriptor, ipt->descriptor+64);
	}
}*/