/** @file    random.h
 ** @brief   Random number generator
 ** @author  Andrea Vedaldi
 **/

/* AUTORIGHTS
Copyright (C) 2007-10 Andrea Vedaldi and Brian Fulkerson

This file is part of VLFeat, available under the terms of the
GNU GPLv2, or (at your option) any later version.
*/

#ifndef __VL_RANDOM_H__
#define __VL_RANDOM_H__

#include "types.h"
#include <stdio.h>
#include <string.h>

/* Period parameters */
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

/** @brief Random numbber generator state */
typedef struct _VlRand {
  vl_uint32 mt [624] ;
  vl_size mti ;
} VlRand ;

/** @name Setting and reading the state */

void vl_rand_init (VlRand * self) ;
void vl_rand_seed (VlRand * self, vl_uint32 s) ;
void vl_rand_seed_by_array (VlRand * self,
                                      vl_uint32 const key [],
                                      vl_size keySize) ;


/** @name Generate random numbers */

inline vl_uint64 vl_rand_uint64 (VlRand * self) ;
inline vl_int64  vl_rand_int63  (VlRand * self) ;
vl_uint32 vl_rand_uint32 (VlRand * self) ;
inline vl_int32  vl_rand_int31  (VlRand * self) ;
inline double    vl_rand_real1  (VlRand * self) ;
inline double    vl_rand_real2  (VlRand * self) ;
inline double    vl_rand_real3  (VlRand * self) ;
inline double    vl_rand_res53  (VlRand * self) ;
inline vl_uindex vl_rand_uindex (VlRand * self, vl_uindex range) ;

/* VL_RANDOM_H */
#endif
