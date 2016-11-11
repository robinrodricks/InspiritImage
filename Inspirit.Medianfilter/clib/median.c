#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "AS3.h"
 
void ctmf( 
			const unsigned char *src, unsigned char *dst, 
			int width, int height, int src_step_row, int dst_step_row, 
			int r, int channels, unsigned long memsize );
void ctmfAsync ( 
			const unsigned char *const src, unsigned char *const dst,
			const int width, const int height,
			const int src_step, const int dst_step,
			const int r, const int cn, const long unsigned int memsize, int *progress );

#include "ctmf.c"

		


unsigned char *result;
unsigned char *source;

int width = 320;
int height = 240;
int area;
int radius = 3;
int channels = 3;

int asyncProgress;

static AS3_Val setupMedianFilter(void* self, AS3_Val args)
{
	AS3_ArrayValue(args, "IntType, IntType, IntType, IntType", &width, &height, &channels, &radius );
	
	area = width * height * (channels == 3 ? 4 : channels);
	
	if(result)	free(result);
	if(source)	free(source);
	
	result = (unsigned char*)malloc( area );
	source = (unsigned char*)malloc( area );
	
	return 0;
}


static AS3_Val runMedianFilter(void* self, AS3_Val args)
{
	AS3_ArrayValue(args, "IntType", &radius );
	
	ctmf ( source, result, width, height, width*channels, width*channels, radius, channels,	256 * 1024 );
	
	return 0;
}

static AS3_Val runMedianFilterAsync(void* self, AS3_Val args)
{
	AS3_ArrayValue(args, "IntType", &radius );
	
	asyncProgress = 0;
	
	ctmfAsync ( source, result, width, height, width*channels, width*channels, radius, channels, 64 * 1024, &asyncProgress );
	
	return 0;
}

static AS3_Val dispose(void* self, AS3_Val args)
{
	free(result);
	free(source);
	
	return 0;
}

static AS3_Val getPointers(void* self, AS3_Val args)
{
	AS3_Val pointers = AS3_Array("AS3ValType", NULL);
	
	AS3_Set(pointers, AS3_Int(0), AS3_Ptr(&result));
	AS3_Set(pointers, AS3_Int(1), AS3_Ptr(&source));
	AS3_Set(pointers, AS3_Int(2), AS3_Ptr(&asyncProgress));
	
	return pointers;
}

int main()
{
	AS3_Val setupMedianFilter_m = AS3_Function( NULL, setupMedianFilter );
	AS3_Val runMedianFilter_m = AS3_Function( NULL, runMedianFilter );
	AS3_Val dispose_m = AS3_Function( NULL, dispose );
	AS3_Val runMedianFilterAsync_m = AS3_FunctionAsync( NULL, runMedianFilterAsync );
	AS3_Val getPointers_m = AS3_Function( NULL, getPointers );

	AS3_Val result = AS3_Object("setupMedianFilter: AS3ValType, runMedianFilter: AS3ValType, dispose: AS3ValType, runMedianFilterAsync: AS3ValType, getPointers: AS3ValType",
								setupMedianFilter_m, runMedianFilter_m, dispose_m, runMedianFilterAsync_m, getPointers_m );

	AS3_Release( setupMedianFilter_m );
	AS3_Release( runMedianFilter_m );
	AS3_Release( dispose_m );
	AS3_Release( runMedianFilterAsync_m );
	AS3_Release( getPointers_m );

	AS3_LibInit( result );

	return 0;
}