void homography_from_4pt(const double x1, const double x2, const double y1, const double y2, 
						const double z1, const double z2, const double w1, const double w2, double *cgret);
void homography_from_4corresp(double *R,
									const double u1, const double v1, const double up1, const double vp1,
                                    const double u2, const double v2, const double up2, const double vp2,
                                    const double u3, const double v3, const double up3, const double vp3,
                                    const double u4, const double v4, const double up4, const double vp4);

void set_bottom_right_coefficient_to_one(double *H);
void normalizePoints(const int number_of_correspondences, const double *u_v_up_vp, double *normalized_u_v_up_vp, double *T1, double *T2inv);


void homography_from_4pt(const double x1, const double x2, const double y1, const double y2, 
						const double z1, const double z2, const double w1, const double w2, double *cgret)
{
	double t1 = x1;
	double t2 = z1;
	double t4 = y2;
	double t5 = t1 * t2 * t4;
	double t6 = w2;
	double t7 = t1 * t6;
	double t8 = t2 * t7;
	double t9 = z2;
	double t10 = t1 * t9;
	double t11 = y1;
	double t14 = x2;
	double t15 = w1;
	double t16 = t14 * t15;
	double t18 = t16 * t11;
	double t20 = t15 * t11 * t9;
	double t21 = t15 * t4;
	double t24 = t15 * t9;
	double t25 = t2 * t4;
	double t26 = t6 * t2;
	double t27 = t6 * t11;
	double t28 = t9 * t11;
	double t30 = (double)0.1e1 / (-t24 + t21 - t25 + t26 - t27 + t28);
	double t32 = t1 * t15;
	double t35 = t14 * t11;
	double t41 = t4 * t1;
	double t42 = t6 * t41;
	double t43 = t14 * t2;
	double t46 = t16 * t9;
	double t48 = t14 * t9 * t11;
	double t51 = t4 * t6 * t2;
	double t55 = t6 * t14;
	cgret[0] = -(-t5 + t8 + t10 * t11 - t11 * t7 - t16 * t2 + t18 - t20 + t21 * t2) * t30;
	cgret[1] = (t5 - t8 - t32 * t4 + t32 * t9 + t18 - t2 * t35 + t27 * t2 - t20) * t30;
	cgret[2] = t1;
	cgret[3] = (-t9 * t7 + t42 + t43 * t4 - t16 * t4 + t46 - t48 + t27 * t9 - t51) * t30;
	cgret[4] = (-t42 + t41 * t9 - t55 * t2 + t46 - t48 + t55 * t11 + t51 - t21 * t9) * t30;
	cgret[5] = t14;
	cgret[6] = (-t10 + t41 + t43 - t35 + t24 - t21 - t26 + t27) * t30;
	cgret[7] = (-t7 + t10 + t16 - t43 + t27 - t28 - t21 + t25) * t30;
	//cgret[8] = 1;
}

void homography_from_4corresp(double *R,
									const double u1, const double v1, const double up1, const double vp1,
                                    const double u2, const double v2, const double up2, const double vp2,
                                    const double u3, const double v3, const double up3, const double vp3,
                                    const double u4, const double v4, const double up4, const double vp4)
{
	double Hr[9], Hl[9];

	homography_from_4pt(u1, v1, u2, v2, u3, v3, u4, v4, &*Hr);
	homography_from_4pt(up1, vp1, up2, vp2, up3, vp3, up4, vp4, &*Hl);

	// the following code computes R = Hl * inverse Hr
	double t2 = Hr[4]-Hr[7]*Hr[5];
	double t4 = Hr[0]*Hr[4];
	double t5 = Hr[0]*Hr[5];
	double t7 = Hr[3]*Hr[1];
	double t8 = Hr[2]*Hr[3];
	double t10 = Hr[1]*Hr[6];
	double t12 = Hr[2]*Hr[6];
	double t15 = (double)1.0 / (t4-t5*Hr[7]-t7+t8*Hr[7]+t10*Hr[5]-t12*Hr[4]);
	double t18 = -Hr[3]+Hr[5]*Hr[6];
	double t23 = -Hr[3]*Hr[7]+Hr[4]*Hr[6];
	double t28 = -Hr[1]+Hr[2]*Hr[7];
	double t31 = Hr[0]-t12;
	double t35 = Hr[0]*Hr[7]-t10;
	double t41 = -Hr[1]*Hr[5]+Hr[2]*Hr[4];
	double t44 = t5-t8;
	double t47 = t4-t7;
	double t48 = t2*t15;
	double t49 = t28*t15;
	double t50 = t41*t15;
	R[0] = Hl[0]*t48+Hl[1]*(t18*t15)-Hl[2]*(t23*t15);
	R[1] = Hl[0]*t49+Hl[1]*(t31*t15)-Hl[2]*(t35*t15);
	R[2] = -Hl[0]*t50-Hl[1]*(t44*t15)+Hl[2]*(t47*t15);
	R[3] = Hl[3]*t48+Hl[4]*(t18*t15)-Hl[5]*(t23*t15);
	R[4] = Hl[3]*t49+Hl[4]*(t31*t15)-Hl[5]*(t35*t15);
	R[5] = -Hl[3]*t50-Hl[4]*(t44*t15)+Hl[5]*(t47*t15);
	R[6] = Hl[6]*t48+Hl[7]*(t18*t15)-t23*t15;
	R[7] = Hl[6]*t49+Hl[7]*(t31*t15)-t35*t15;
	R[8] = -Hl[6]*t50-Hl[7]*(t44*t15)+t47*t15;
}

void set_bottom_right_coefficient_to_one(double *H)
{
	const double inv_H33 = (double)1.0 / H[8];
	H[0] *= inv_H33;
	H[1] *= inv_H33;
	H[2] *= inv_H33;
	H[3] *= inv_H33;
	H[4] *= inv_H33;
	H[5] *= inv_H33;
	H[6] *= inv_H33;
	H[7] *= inv_H33;
	H[8] = (double)1.0;
}

void normalizePoints(const int number_of_correspondences, const double *u_v_up_vp, double *normalized_u_v_up_vp, double *T1, double *T2inv)
{
	const double invN = (double)1.0 / (double)number_of_correspondences;
	const double sqrt2N = SQRT2 * (double)number_of_correspondences;
	
	double u_sum = 0., v_sum = 0., up_sum = 0., vp_sum = 0;
	int i, j;
	
	for(i = 0, j = 0; i < number_of_correspondences; i++) {
		u_sum  += u_v_up_vp[j++];
		v_sum  += u_v_up_vp[j++];
		up_sum += u_v_up_vp[j++];
		vp_sum += u_v_up_vp[j++];
	}

	double u_mean  = u_sum  * invN;
	double v_mean  = v_sum  * invN;
	double up_mean = up_sum * invN;
	double vp_mean = vp_sum * invN;

	// translate mean to origin, compute sum of distances from origin
	double dist_sum = 0, distp_sum = 0;
	double n_up_vp_um, n_up_vp_vm;
	for(i = 0; i < number_of_correspondences; i++) 
	{
		normalized_u_v_up_vp[4 * i    ] = ( n_up_vp_um = u_v_up_vp[4 * i    ] - u_mean );
		normalized_u_v_up_vp[4 * i + 1] = ( n_up_vp_vm = u_v_up_vp[4 * i + 1] - v_mean );

		dist_sum += fast_sqrt( dSquare(n_up_vp_um) + dSquare(n_up_vp_vm) );

		normalized_u_v_up_vp[4 * i + 2] = ( n_up_vp_um = u_v_up_vp[4 * i + 2] - up_mean );
		normalized_u_v_up_vp[4 * i + 3] = ( n_up_vp_vm = u_v_up_vp[4 * i + 3] - vp_mean );

		distp_sum += fast_sqrt( dSquare(n_up_vp_um) + dSquare(n_up_vp_vm) );
	}

	// compute normalizing scale factor ( average distance from origin = sqrt(2) )
	double scale  = sqrt2N / dist_sum;
	double scalep = sqrt2N / distp_sum;

	// apply scaling
	for(i = 0, j = 0; i < number_of_correspondences; i++) {
		normalized_u_v_up_vp[j++] *= scale;
		normalized_u_v_up_vp[j++] *= scale;
		normalized_u_v_up_vp[j++] *= scalep;
		normalized_u_v_up_vp[j++] *= scalep;
	}

	// assemble transformation Matrices, used at denormalization

	T1[1] = T1[3] = T1[6] = T1[7] = 0.0;
	T1[0] = scale;
	T1[2] = -scale * u_mean;
	T1[4] = scale;
	T1[5] = -scale * v_mean;
	T1[8] = (double)1.0;
	
	T2inv[1] = T2inv[3] = T2inv[6] = T2inv[7] = 0.0;
	T2inv[0] = (double)1. / scalep;
	T2inv[2] = up_mean;
	T2inv[4] = (double)1. / scalep;
	T2inv[5] = vp_mean;
	T2inv[8] = (double)1.0;
}

void denormalizeHomography(double *T1, double *T2inv, double *H)
{
	double tmp[9];
	MultiplyMat(T2inv, H, &*tmp, 3, 3, 3, 3);
	MultiplyMat(tmp, T1, H, 3, 3, 3, 3);
}

void get_4_prosac_indices(const int n_max, int *n1, int *n2, int *n3, int *n4)
{
	*n1 = rand() % n_max;
	do *n2 = rand() % n_max; while(*n2 == *n1);
	do *n3 = rand() % n_max; while(*n3 == *n1 || *n3 == *n2);
	do *n4 = rand() % n_max; while(*n4 == *n1 || *n4 == *n2 || *n4 == *n3);
	
	*n1 <<= 2;
	*n2 <<= 2;
	*n3 <<= 2;
	*n4 <<= 2;
}

int compute_inliers(const IPointMatch *matches, const int number_of_correspondences, const double *H, int *inlier_ids, const double threshold)
{
	int i, num_inliers = 0;
	double x, y, Z;
	double eup, evp;
	
	const double m11 = H[0], m12 = H[1], m13 = H[2], m21 = H[3], m22 = H[4], m23 = H[5], m31 = H[6], m32 = H[7], m33 = H[8];

	for(i = 0; i < number_of_correspondences; i++) 
	{
		x = (double)matches[i].second->x, y = (double)matches[i].second->y;
		Z = (double)1.0 / (m31*x + m32*y + m33);
		
		eup = (m11*x + m12*y + m13) * Z;
		evp = (m21*x + m22*y + m23) * Z;
		
		double distance = dSquare((double)matches[i].first->x - eup) + dSquare((double)matches[i].first->y - evp);
		if(distance < threshold)
		{
			inlier_ids[num_inliers] = i;
			num_inliers++;
		}
  }

  return num_inliers;
}

int nice_homography(const double * H)
{
  double det = H[0] * H[4] - H[3] * H[1];
  if (det < 0) return 0;
  double N1 = fast_sqrt( dSquare(H[0]) + dSquare(H[3]) );
  if (N1 > 4) return 0;
  if (N1 < 0.1) return 0;
  double N2 = fast_sqrt( dSquare(H[1]) + dSquare(H[4]) );
  if (N2 > 4) return 0;
  if (N2 < 0.1) return 0;
  double N3 = fast_sqrt( dSquare(H[6]) + dSquare(H[7]) );
  if (N3 > 0.002) return 0;
  return 1;
}