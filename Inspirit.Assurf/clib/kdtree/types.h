#ifndef __KD_TYPES__
#define __KD_TYPES__


#define VL_NAN_F 0x7FC00000UL
#define VL_INFINITY_F 0x7F800000UL
#define VL_NAN_D 0x7FF8000000000000ULL
#define VL_INFINITY_D 0x7FF0000000000000ULL

#define VL_MIN(x,y) (((x)<(y))?(x):(y))
#define VL_MAX(x,y) (((x)>(y))?(x):(y))

#define VL_CAT(x,y) x ## y
#define VL_XCAT(x,y) VL_CAT(x,y)
#define VL_XCAT3(x,y,z) VL_XCAT(VL_XCAT(x,y),z)
#define VL_XCAT4(x,y,z,u) VL_XCAT(VL_XCAT3(x,y,z),u)
#define VL_XCAT5(x,y,z,u,v) VL_XCAT(VL_XCAT4(x,y,z,u),v)

typedef long long           vl_int64 ;
typedef int                 vl_int32 ;
typedef short               vl_int16 ;
typedef char                vl_int8  ;
typedef long long unsigned  vl_uint64 ;
typedef int       unsigned  vl_uint32 ;
typedef short     unsigned  vl_uint16 ;
typedef char      unsigned  vl_uint8 ;
typedef int                 vl_int ;
typedef unsigned int        vl_uint ;
typedef int                 vl_bool ;
typedef vl_int32            vl_intptr ;
typedef vl_uint32           vl_uintptr ;
typedef vl_uint32           vl_size ;
typedef vl_int32            vl_index ;
typedef vl_uint32           vl_uindex ;

#define VL_EXPORT

#define VL_TYPE_FLOAT   1
#define VL_TYPE_DOUBLE  2
#define VL_TYPE_INT8    3
#define VL_TYPE_UINT8   4
#define VL_TYPE_INT16   5
#define VL_TYPE_UINT16  6
#define VL_TYPE_INT32   7
#define VL_TYPE_UINT32  8
#define VL_TYPE_INT64   9
#define VL_TYPE_UINT64  10

typedef vl_uint32 vl_type ;

#endif
