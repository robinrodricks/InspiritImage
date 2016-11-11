#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "AS3.h"

#include "utils.h"
#include "homography.c"

unsigned char *pyr_img1;
unsigned char *pyr_img2;
unsigned char *pyr_img3;
unsigned char *pyr_img_blur1;
unsigned char *pyr_img_blur2;
unsigned char *pyr_img_blur3;
unsigned char *pyr_img_prev1;
unsigned char *pyr_img_prev2;
unsigned char *pyr_img_prev3;
unsigned char *img_mask;

int pyr_img1_width = 320;
int pyr_img1_height = 240;
int pyr_img2_width = 320;
int pyr_img2_height = 240;
int pyr_img3_width = 320;
int pyr_img3_height = 240;

unsigned char *pyr_img_prev[3];
unsigned char *pyr_img[3];
unsigned char *pyr_img_blur[3];

double *descriptorsRef;
double *descriptorsCurr;
double *arcamera;

int *result_b;
int *currRefIndexes;

int max_ref_points_pool = 5000;
int max_screen_points = 2000;
int max_points_pool = 500;
int max_objects_pool = 100;

IPoint *referencePoints;
IPoint *screenPoints;
IPoint *screenPointsPrev;
IPointMatch *frameMatches;
IPointMatch *refMatches;
RefObject *refObjectsMap;

int referencePointsCount = 0;
int referenceCount = 0;
int screenPointsCount = 0;
int prevFramePointsCount = 0;
int matchedPointsCount = 0;
int foundReferenceCount = 0;

int screenPointsThresh = 10;
double screenPointsShitomasi = 70.0;
int useMask = 0;
int supressNeighbors = 0;
double supressDist = 15.0;

int DESCRIPTOR_SIZE = 36;
int PATCH_SIZE = 64;
int DETECT_PRECISION = 2;


// --------------------------------------------
// KDTREE STRUCTURE

#include "match_kdtree.c"

// --------------------------------------------
// Points Detectors and Descriptors

#include "detector_fast.c"
#include "descriptor_surf.c"
#include "tracker.c"

// --------------------------------------------
// 3D POSE ESTIMATION

double estimatePoseFromMatches(double *camera_info, IPointMatch *matches, double model_matrix[12], const int matchesCount, const int width, const int height);
double estimatePoseFromCorners(double *camera_info, RefObject *object);
double estimatePoseFromPoints(double *camera_info, double *scr_pts, double *ref_pts, double model_matrix[12], const int count, const int width, const int height);

// --------------------------------------------

// --------------------------------------------
// EXPORT POINTS DATA

#include "export.c"

// --------------------------------------------

/*static void * custom_memmove( void * destination, const void * source, size_t num ) {
 
  void *result;
  __asm__("%0 memmove(%1, %2, %3)\n" : "=r"(result) : "r"(destination), "r"(source), "r"(num));
  return result;
}*/
 
 
static void * custom_memcpy ( void * destination, const void * source, size_t num ) {
  void *result;
 
  __asm__("%0 memcpy(%1, %2, %3)\n" : "=r"(result) : "r"(destination), "r"(source), "r"(num));
  return result;
}
 
 
 
static void * custom_memset ( void * ptr, int value, size_t num ) {
  void *result;
  __asm__("%0 memset(%1, %2, %3)\n" : "=r"(result) : "r"(ptr), "r"(value), "r"(num));
  return result;
}
 
 
//#define memmove custom_memmove
#define memcpy custom_memcpy
#define memset custom_memset

// --------------------------------------------

AS3_Val setupGlobalBuffers(void* self, AS3_Val args)
{	
	AS3_ArrayValue(args, "IntType, IntType, IntType, IntType", &max_ref_points_pool, &max_objects_pool, &max_points_pool, &DETECT_PRECISION );
	
	if(DETECT_PRECISION == 0) DESCRIPTOR_SIZE = 128;
	if(DETECT_PRECISION == 1) DESCRIPTOR_SIZE = 64;
	if(DETECT_PRECISION == 2) DESCRIPTOR_SIZE = 36;
	
	if(referencePoints) free(referencePoints);
	if(frameMatches) free(frameMatches);
	if(refMatches) free(refMatches);
	if(refObjectsMap) free(refObjectsMap);
	if(descriptorsRef) free(descriptorsRef);
	if(arcamera) free(arcamera);
	if(result_b) free(result_b);
	if(currRefIndexes) free(currRefIndexes);
	
	referencePoints = (IPoint*)malloc( (max_ref_points_pool + max_screen_points + max_points_pool*2) * sizeof(IPoint) );
	screenPoints = referencePoints+max_ref_points_pool;
	screenPointsPrev = screenPoints+max_screen_points;
	frameMatches = (IPointMatch*)malloc( max_points_pool*2 * sizeof(IPointMatch) );
	refMatches = (IPointMatch*)malloc( max_objects_pool*max_points_pool * sizeof(IPointMatch) );
	
	refObjectsMap = (RefObject *)malloc( max_objects_pool * sizeof(RefObject) );
	
	descriptorsRef = (double*)malloc( ((DESCRIPTOR_SIZE * (max_ref_points_pool+max_points_pool)) ) * sizeof(double) );
	descriptorsCurr = descriptorsRef + (DESCRIPTOR_SIZE * max_ref_points_pool);
	
	arcamera = (double*)malloc(4 * sizeof(double));
	
	memset(refObjectsMap, 0, max_objects_pool * sizeof(refObjectsMap));
	
	result_b = (int *)malloc((max_points_pool + 10) * sizeof(int));
	currRefIndexes = (int *)malloc(max_objects_pool * sizeof(int));
	
	calcGaussWeights( (DETECT_PRECISION == 2) ? 15 : 20 );

	return AS3_Ptr(result_b);
}

AS3_Val setupImageHolders(void* self, AS3_Val args)
{	
	AS3_ArrayValue(args, "IntType, IntType", &pyr_img1_width, &pyr_img1_height );
	
	if(pyr_img1) free(pyr_img1);
	
	int off = pyr_img1_width * pyr_img1_height;
	
	pyr_img1 = (unsigned char*)malloc( off * 2 * sizeof(unsigned char));
	pyr_img_blur1 = pyr_img1+off;
	
	pyr_img[0] = pyr_img1;	
	pyr_img_blur[0] = pyr_img_blur1;
	
	AS3_Val pointers = AS3_Array("AS3ValType", NULL);
	
	AS3_Set(pointers, AS3_Int(0), AS3_Ptr(pyr_img1) );
	AS3_Set(pointers, AS3_Int(1), AS3_Ptr(pyr_img_blur1) );
	
	return pointers;
}

AS3_Val setupImagePyramid(void* self, AS3_Val args)
{	
	AS3_ArrayValue(args, "IntType, IntType", &pyr_img1_width, &pyr_img1_height );
	
	if(pyr_img1) free(pyr_img1);
	
	int off = pyr_img1_width * pyr_img1_height;
	
	pyr_img1 = (unsigned char*)malloc( off * 4 * sizeof(unsigned char));
	pyr_img_blur1 = pyr_img1+off;
	pyr_img_prev1 = pyr_img_blur1+off;
	img_mask = pyr_img_prev1+off;
	
	//
	
	if(pyr_img2) free(pyr_img2);
	
	pyr_img2_width = pyr_img1_width >> 1;
	pyr_img2_height = pyr_img1_height >> 1;
	off = pyr_img2_width * pyr_img2_height;
	
	pyr_img2 = (unsigned char*)malloc( off * 3 * sizeof(unsigned char));
	pyr_img_blur2 = pyr_img2 + off;
	pyr_img_prev2 = pyr_img_blur2 + off;
	
	//
	
	if(pyr_img3) free(pyr_img3);
	
	pyr_img3_width = pyr_img1_width >> 2;
	pyr_img3_height = pyr_img1_height >> 2;
	
	off = pyr_img3_width * pyr_img3_height;
	
	pyr_img3 = (unsigned char*)malloc( off * 3 * sizeof(unsigned char));
	pyr_img_blur3 = pyr_img3 + off;
	pyr_img_prev3 = pyr_img_blur3 + off;
	
	//
	
	pyr_img_prev[0] = pyr_img_prev1;
	pyr_img_prev[1] = pyr_img_prev2;
	pyr_img_prev[2] = pyr_img_prev3;
	
	pyr_img[0] = pyr_img1;
	pyr_img[1] = pyr_img2;
	pyr_img[2] = pyr_img3;
	
	pyr_img_blur[0] = pyr_img_blur1;
	pyr_img_blur[1] = pyr_img_blur2;
	pyr_img_blur[2] = pyr_img_blur3;
	
	int sizeW[3] = {pyr_img1_width, pyr_img2_width, pyr_img3_width};
	int sizeH[3] = {pyr_img1_height, pyr_img2_height, pyr_img3_height};
	initOpticalFlow( pyr_img_prev, pyr_img, sizeW, sizeH, 16, CV_LKFLOW_INITIAL_GUESSES );
	
	AS3_Val pointers = AS3_Array("AS3ValType", NULL);
	
	AS3_Set(pointers, AS3_Int(0), AS3_Ptr(pyr_img1) );
	AS3_Set(pointers, AS3_Int(1), AS3_Ptr(pyr_img2) );
	AS3_Set(pointers, AS3_Int(2), AS3_Ptr(pyr_img3) );
	AS3_Set(pointers, AS3_Int(3), AS3_Ptr(pyr_img_blur1) );
	AS3_Set(pointers, AS3_Int(4), AS3_Ptr(pyr_img_blur2) );
	AS3_Set(pointers, AS3_Int(5), AS3_Ptr(pyr_img_blur3) );
	AS3_Set(pointers, AS3_Int(6), AS3_Ptr(img_mask) );
	
	return pointers;
}

AS3_Val createRefObject(void* self, AS3_Val args)
{
	int iw = 320;
	int ih = 240;
	
	AS3_ArrayValue(args, "IntType, IntType", &iw, &ih );
	
	RefObject obj;
	
	obj.index = referenceCount;
	obj.width = iw;
	obj.height = ih;
	obj.pointsCount = 0;
	obj.descriptors = descriptorsRef+(referencePointsCount*DESCRIPTOR_SIZE);
	obj.points = referencePoints+referencePointsCount;
	obj.matches = refMatches + ( referenceCount * max_points_pool );
	
	refObjectsMap[referenceCount] = obj;
	
	RefObject *obj_ptr = &refObjectsMap[referenceCount];
	
	foundReferenceCount = referenceCount++;
	
	AS3_Val pointers = AS3_Array("AS3ValType", NULL);
	
	AS3_Set(pointers, AS3_Int(0), AS3_Int(obj_ptr->index));
	AS3_Set(pointers, AS3_Int(1), AS3_Ptr(&obj_ptr->pointsCount));
	AS3_Set(pointers, AS3_Int(2), AS3_Ptr(&obj_ptr->matchedPointsCount));
	AS3_Set(pointers, AS3_Int(3), AS3_Ptr(&obj_ptr->poseError));
	AS3_Set(pointers, AS3_Int(4), AS3_Ptr(&*obj_ptr->homography));
	AS3_Set(pointers, AS3_Int(5), AS3_Ptr(&*obj_ptr->pose));
	
	return pointers;
}

AS3_Val pushDataToRefObject(void* self, AS3_Val args)
{
	int refID = 0;
	int refObjPointsCount = 0;
	AS3_Val src;
	FILE *input;
	int *iptr;
	double *dptr;
	
	AS3_ArrayValue( args, "IntType, IntType, AS3ValType", &refID, &refObjPointsCount, &src );
	
	RefObject *obj = &refObjectsMap[refID];
	
	input = funopen((void *)src, readba, writeba, seekba, closeba);
	
	int chunck_size = ( (DESCRIPTOR_SIZE * 8) + 40 );
	char out[chunck_size];
	
	int i, j;
	for(i = 0; i < refObjPointsCount; i++)
	{
		fread(out, 1, chunck_size, input);
		iptr = (int*)out;
		
		IPoint *ipt = obj->points + i;
			
		ipt->index = referencePointsCount;
		iptr++;
		ipt->refIndex = refID;
		iptr++;
		ipt->pos = *iptr++;
		ipt->x = *iptr++;
		ipt->y = *iptr++;
		/*ipt->scale = */*iptr++;
		ipt->localIndex = i;
		
		dptr = (double*)(out+24);
			
		ipt->score = *dptr++;
		ipt->orientation = *dptr++;
		
		ipt->descriptor = descriptorsRef+(referencePointsCount*DESCRIPTOR_SIZE);
		
		j = DESCRIPTOR_SIZE;
		double *descr = ipt->descriptor;
		while( --j > -1 )
		{
			*descr++ = *dptr++;
		}
		
		referencePointsCount ++;
	}
	
	obj->pointsCount = refObjPointsCount;
	
	fclose( input );
	
	return 0;
}

AS3_Val pushImageToRefObject(void* self, AS3_Val args)
{
	int count;
	
	int refID = 0;
	int iw = 320;
	int ih = 240;
	int maxPoints = 300;
	double scale = 1.0;
	int i, j;
	IPoint *pts;
	
	AS3_ArrayValue(args, "IntType, IntType, IntType, IntType, DoubleType", &refID, &iw, &ih, &maxPoints, &scale );
	
	int ip_thresh = ((supressNeighbors == 0) ? 45 : 15);
	
	count = detectPointsLevel(0, iw, ih, referencePoints+referencePointsCount, ip_thresh, 70.0, maxPoints, 25, max_ref_points_pool - referencePointsCount);
	
	if(supressNeighbors) 
	{
		supressDist = (double)12.0 / scale;
		count = supressNeighborPoints(referencePoints+referencePointsCount, count, supressDist * supressDist);
	}
	
	if(DESCRIPTOR_SIZE == 36)
	{
		fastSURFDescriptors( referencePoints+referencePointsCount, count, descriptorsRef+(referencePointsCount*DESCRIPTOR_SIZE) );
	}
	else if(DESCRIPTOR_SIZE == 64)
	{
		SURFDescriptors( referencePoints+referencePointsCount, count, descriptorsRef+(referencePointsCount*DESCRIPTOR_SIZE), 0 );
	}
	else if(DESCRIPTOR_SIZE == 128)
	{
		SURFDescriptors( referencePoints+referencePointsCount, count, descriptorsRef+(referencePointsCount*DESCRIPTOR_SIZE), 1 );
	}
	
	pts = referencePoints+referencePointsCount;
	
	for(i = 0, j = referencePointsCount; i < count; i++, j++)
	{
		pts[i].x *= scale;
		pts[i].y *= scale;
		pts[i].refIndex = refID;
		pts[i].index = j;
		pts[i].localIndex = i;
	}
	
	RefObject *obj_ptr = &refObjectsMap[refID];
	obj_ptr->pointsCount += count;	
	
	referencePointsCount += count;
	
	return 0;
}

AS3_Val clearReferenceObjects(void* self, AS3_Val args)
{
	int i;
	for(i = 0; i < referenceCount; i++)
	{
		if(refObjectsMap[i].kdf) vl_kdforest_delete(refObjectsMap[i].kdf);
	}
	memset(refObjectsMap, 0, max_objects_pool * sizeof(RefObject));
	
	referencePointsCount = 0;
	referenceCount = 0;
	
	return 0;
}

AS3_Val buildRefIndex(void* self, AS3_Val args)
{		
	int i;
	for(i=0; i<referenceCount; i++)
	{
		RefObject *obj = &refObjectsMap[i];
		obj->kdf = vl_kdforest_new ( DESCRIPTOR_SIZE, 5 );
		vl_kdforest_build (obj->kdf, obj->pointsCount, obj->descriptors);
		vl_kdforest_set_max_num_comparisons(obj->kdf, 81);
	}
	
	result_b[max_points_pool+3] = referencePointsCount;
	supressNeighbors = 0;
	
	return 0;
}

void relocateFeatures(IPointMatch *matches, const int count)
{	
	int i;	
	for(i = 0; i < count; i++)
	{		
		memcpy(screenPointsPrev + i, matches[i].first, sizeof(IPoint));
		matches[i].first = &screenPointsPrev[i];
		matches[i].first->matched = 0;
	}
	
	prevFramePointsCount = count;
}

AS3_Val runTask(void* self, AS3_Val args)
{	
	int maxPoints = 300;
	int objID = 0;
	int options = 0;
	int count = 0;
	int skiped = 0;
	
	AS3_ArrayValue(args, "IntType, IntType, IntType, IntType", &maxPoints, &useMask, &objID, &options );
	
	if(screenPointsCount > maxPoints - 50)
	{
		if(screenPointsThresh < 150) screenPointsThresh += 1;
		if(screenPointsThresh >= 150 && screenPointsShitomasi < 150) screenPointsShitomasi += 1;
	}
	else if(screenPointsCount < 100)
	{
		if(screenPointsShitomasi <= 50 && screenPointsThresh > 10) screenPointsThresh -= 1;
		if(screenPointsShitomasi > 50) screenPointsShitomasi -= 1;
	}
	
	screenPointsCount = detectPointsLevel2(3, screenPoints, screenPointsThresh, screenPointsShitomasi, maxPoints, 15, max_screen_points);
	
	//prevFramePointsCount = findPrevFrameMatches(screenPoints, screenPointsCount, frameMatches, prevFramePointsCount);
	prevFramePointsCount = NCCLKTrack(screenPoints, screenPointsCount, frameMatches, prevFramePointsCount);
	
	int goodMatches = 0;
	int i;
	for(i = 0; i < prevFramePointsCount; i++)
	{
		if(frameMatches[i].wasGood == 1) goodMatches++;
	}
	
	RefObject *obj = &refObjectsMap[objID];
	
	IPointMatch origMatches[maxPoints + prevFramePointsCount];
	int origMatCount = 0;
	
	if(goodMatches < 15)
	{
	
		calculateDescriptors:
	
		if(DESCRIPTOR_SIZE == 36)
		{
			fastSURFDescriptors( screenPoints, screenPointsCount, descriptorsCurr );
		}
		else if(DESCRIPTOR_SIZE == 64)
		{
			SURFDescriptors( screenPoints, screenPointsCount, descriptorsCurr, 0 );
		}
		else if(DESCRIPTOR_SIZE == 128)
		{
			SURFDescriptors( screenPoints, screenPointsCount, descriptorsCurr, 1 );
		}

		count = getMatchesByTree(obj, screenPoints, referencePoints, screenPointsCount, referencePointsCount, frameMatches, prevFramePointsCount, 0.7);
		/*count = getMatchesByTree(obj, screenPoints, referencePoints, screenPointsCount, referencePointsCount, frameMatches+prevFramePointsCount, 0, 0.7);
		memcpy(&*origMatches, frameMatches+prevFramePointsCount, count * sizeof(IPointMatch));
		
		prevFramePointsCount = NCCLKTrack(screenPoints, screenPointsCount, frameMatches, prevFramePointsCount);
		memcpy(&*(origMatches+count), frameMatches, prevFramePointsCount * sizeof(IPointMatch));
		memcpy(&*(frameMatches+prevFramePointsCount), &*origMatches, count * sizeof(IPointMatch));
		count += prevFramePointsCount;*/
		
		origMatCount = count;
		memcpy(&*origMatches, frameMatches, count * sizeof(IPointMatch));

		result_b[max_points_pool] = screenPointsCount;
		result_b[max_points_pool+1] = count;
		result_b[max_points_pool+2] = prevFramePointsCount;
	
		count = filterOutliersByAngle(frameMatches, count);
		//count = filterOutliersByLines(frameMatches, count);
		
		locatePlanarObject(frameMatches, &count, &*obj->homography);
	}
	else
	{
		result_b[max_points_pool] = screenPointsCount;
		result_b[max_points_pool+1] = prevFramePointsCount;
		result_b[max_points_pool+2] = prevFramePointsCount;
		count = prevFramePointsCount;
		
		origMatCount = count;
		memcpy(&*origMatches, frameMatches, count * sizeof(IPointMatch));
		
		locatePlanarObject(frameMatches, &count, &*obj->homography);
		
		skiped = 1;
	}
	
	if(skiped == 1 && count < 8)
	{
		skiped = 0;
		//prevFramePointsCount = count;
		memcpy(frameMatches, &*origMatches, origMatCount * sizeof(IPointMatch));
		goto calculateDescriptors;
	}
	
	obj->matchedPointsCount = count;
	
	if(options == 2 && count > 4)
	{
		obj->poseError = estimatePoseFromMatches(arcamera, frameMatches, &*obj->pose, count, obj->width, obj->height);
		
		if(skiped == 1 && obj->poseError > 4)
		{
			skiped = 0;
			//prevFramePointsCount = count;
			memcpy(frameMatches, &*origMatches, origMatCount * sizeof(IPointMatch));
			goto calculateDescriptors;
		}
	}	

	memcpy(pyr_img_prev1, pyr_img_blur1, pyr_img1_width * pyr_img1_height);
	memcpy(pyr_img_prev2, pyr_img_blur2, pyr_img2_width * pyr_img2_height);
	memcpy(pyr_img_prev3, pyr_img_blur3, pyr_img3_width * pyr_img3_height);
	//memcpy(pyr_img_prev1, pyr_img1, pyr_img1_width * pyr_img1_height);
	//memcpy(pyr_img_prev2, pyr_img2, pyr_img2_width * pyr_img2_height);
	//memcpy(pyr_img_prev3, pyr_img3, pyr_img3_width * pyr_img3_height);
	
	result_b[max_points_pool+4] = count;
	result_b[max_points_pool+5] = skiped;
	
	if(origMatCount)
	{
		int j;
		for(i = 0; i < origMatCount; i++)
		{
			for(j = 0; j < count; j++)
			{
				if(origMatches[i].second->localIndex == frameMatches[j].second->localIndex){
					origMatches[i].wasGood = frameMatches[j].wasGood;
					break;
				}
			}
		}
		memcpy(frameMatches, &*origMatches, origMatCount * sizeof(IPointMatch));
		relocateFeatures(frameMatches, origMatCount);
	} else {
		relocateFeatures(frameMatches, count);
	}
	
	return AS3_Int(count);
}

AS3_Val runTask2(void* self, AS3_Val args)
{	
	int maxPoints = 300;
	int objNum = 0;
	int options = 0;
	int i, k, count, skiped = 0;
	RefObject *obj;
	
	AS3_ArrayValue(args, "IntType, IntType, IntType", &maxPoints, &useMask, &options );
	
	if(screenPointsCount > maxPoints - 50 /*&& screenPointsThresh < 150*/)
	{
		if(screenPointsThresh < 150) screenPointsThresh += 1;
		if(screenPointsThresh >= 150 && screenPointsShitomasi < 150) screenPointsShitomasi += 1;
	}
	else if(screenPointsCount < 100 /*&& screenPointsThresh > 10*/)
	{
		if(screenPointsShitomasi <= 50 && screenPointsThresh > 10) screenPointsThresh -= 1;
		if(screenPointsShitomasi > 50) screenPointsShitomasi -= 1;
	}
	
	screenPointsCount = detectPointsLevel2(3, screenPoints, screenPointsThresh, screenPointsShitomasi, maxPoints, 15, max_screen_points);
	//prevFramePointsCount = findPrevFrameMatches(screenPoints, screenPointsCount, frameMatches, prevFramePointsCount);
	prevFramePointsCount = NCCLKTrack(screenPoints, screenPointsCount, frameMatches, prevFramePointsCount);
	
	int goodMatches = 0;
	for(i = 0; i < prevFramePointsCount; i++)
	{
		if(frameMatches[i].wasGood == 1) goodMatches++;
	}
	
	if(goodMatches < 15*foundReferenceCount)
	{
	
		calculateDescriptors:
	
		if(DESCRIPTOR_SIZE == 36)
		{
			fastSURFDescriptors( screenPoints, screenPointsCount, descriptorsCurr );
		}
		else if(DESCRIPTOR_SIZE == 64)
		{
			SURFDescriptors( screenPoints, screenPointsCount, descriptorsCurr, 0 );
		}
		else if(DESCRIPTOR_SIZE == 128)
		{
			SURFDescriptors( screenPoints, screenPointsCount, descriptorsCurr, 1 );
		}
		
		count = getMatchesByMultiTree(screenPoints, referencePoints, screenPointsCount, referencePointsCount, frameMatches, prevFramePointsCount, 0.7);
		
		k = 0;
		for(i = 0; i < referenceCount; i++)
		{
			obj = &refObjectsMap[ i ];
			memcpy(frameMatches + k, obj->matches, obj->matchedPointsCount * sizeof(IPointMatch));
			k += obj->matchedPointsCount;
		}
	} 
	else 
	{
		count = prevFramePointsCount;
		skiped = 1;
		sortMatchesByObjects(refObjectsMap, frameMatches, referenceCount, count);
	}
	
	result_b[max_points_pool] = screenPointsCount;
	result_b[max_points_pool+1] = count;
	result_b[max_points_pool+2] = prevFramePointsCount;
	
	foundReferenceCount = 0;
	
	k = 0;
	IPointMatch objMatches[maxPoints];
	for(i = 0; i < referenceCount; i++)
	{
	
		obj = &refObjectsMap[ i ];
		
		if(obj->matchedPointsCount < 5) continue;
		
		currRefIndexes[objNum++] = i;
		
		obj->matchedPointsCount = filterOutliersByAngle(obj->matches, obj->matchedPointsCount);
		//obj->matchedPointsCount = filterOutliersByLines(obj->matches, obj->matchedPointsCount);
		
		locatePlanarObject(obj->matches, &obj->matchedPointsCount, &*obj->homography);
		
		if(skiped == 1 && obj->prevMatchedPointsCount >= 15 && obj->matchedPointsCount < 8)
		{
			skiped = 0;
			goto calculateDescriptors;
		}
		
		if(options == 2 && obj->matchedPointsCount > 4)
		{
			if(obj->matchedPointsCount <= 16)
			{
				obj->poseError = estimatePoseFromMatches(arcamera, obj->matches, &*obj->pose, obj->matchedPointsCount, obj->width, obj->height);
			}
			else
			{
				IPointMatch perpMatches[8];
				perpendicularRegression(obj->matches, obj->matchedPointsCount, &*perpMatches);
				
				obj->poseError = estimatePoseFromMatches(arcamera, &*perpMatches, &*obj->pose, 8, obj->width, obj->height);
			}
		}
		
		memcpy(&*(objMatches + k), obj->matches, obj->matchedPointsCount * sizeof(IPointMatch));
		k += obj->matchedPointsCount;
		
		if(obj->matchedPointsCount >= 4) foundReferenceCount++;
	}
	
	if(foundReferenceCount == 0) foundReferenceCount = referenceCount;
	
	result_b[max_points_pool+4] = k;
	result_b[max_points_pool+5] = skiped;

	memcpy(pyr_img_prev1, pyr_img_blur1, pyr_img1_width * pyr_img1_height);
	memcpy(pyr_img_prev2, pyr_img_blur2, pyr_img2_width * pyr_img2_height);
	memcpy(pyr_img_prev3, pyr_img_blur3, pyr_img3_width * pyr_img3_height);
	
	int j;
	for(i = 0; i < k; i++)
	{
		for(j = 0; j < count; j++)
		{
			if(objMatches[i].first->index == frameMatches[j].first->index){
				frameMatches[j].wasGood = objMatches[i].wasGood;
				break;
			}
		}
	}
	relocateFeatures(frameMatches, count);
	
	return AS3_Int(objNum);
}

AS3_Val getDataPointers(void* self, AS3_Val args)
{
	AS3_Val pointers = AS3_Array("AS3ValType", NULL);
	
	AS3_Set(pointers, AS3_Int(0), AS3_Ptr(result_b));
	AS3_Set(pointers, AS3_Int(1), AS3_Ptr(pyr_img_blur1));
	AS3_Set(pointers, AS3_Int(2), AS3_Ptr(img_mask));
	AS3_Set(pointers, AS3_Int(3), AS3_Ptr(pyr_img1));
	AS3_Set(pointers, AS3_Int(4), AS3_Ptr(arcamera));
	AS3_Set(pointers, AS3_Int(5), AS3_Ptr(currRefIndexes));
	AS3_Set(pointers, AS3_Int(6), AS3_Ptr(&supressNeighbors));
	
	AS3_Set(pointers, AS3_Int(7), AS3_Ptr(pyr_img2));
	AS3_Set(pointers, AS3_Int(8), AS3_Ptr(pyr_img_blur2));
	AS3_Set(pointers, AS3_Int(9), AS3_Ptr(pyr_img3));
	AS3_Set(pointers, AS3_Int(10), AS3_Ptr(pyr_img_blur3));
	
	return pointers;
}

AS3_Val clearBuffers(void* self, AS3_Val args)
{
	if(pyr_img1) free(pyr_img1);
	if(pyr_img2) free(pyr_img2);
	
	if(descriptorsRef) free(descriptorsRef);
	if(descriptorsCurr) free(descriptorsCurr);
	if(arcamera) free(arcamera);
	if(result_b) free(result_b);
	
	if(referencePoints) free(referencePoints);
	if(screenPoints) free(screenPoints);
	if(screenPointsPrev) free(screenPointsPrev);
	if(frameMatches) free(frameMatches);
	if(refMatches) free(refMatches);
	if(refObjectsMap) free(refObjectsMap);
	if(currRefIndexes) free(currRefIndexes);

	return 0;
}


int main() 
{
	AS3_Val setupGlobalBuffers_ = AS3_Function( NULL, setupGlobalBuffers );
	AS3_Val setupImageHolders_ = AS3_Function( NULL, setupImageHolders );
	AS3_Val createRefObject_ = AS3_Function( NULL, createRefObject );
	AS3_Val pushImageToRefObject_ = AS3_Function( NULL, pushImageToRefObject );
	AS3_Val clearReferenceObjects_ = AS3_Function( NULL, clearReferenceObjects );
	AS3_Val buildRefIndex_ = AS3_Function( NULL, buildRefIndex );
	AS3_Val runTask_ = AS3_Function( NULL, runTask );
	AS3_Val runTask2_ = AS3_Function( NULL, runTask2 );
	AS3_Val getDataPointers_ = AS3_Function( NULL, getDataPointers );
	AS3_Val clearBuffers_ = AS3_Function( NULL, clearBuffers );
	AS3_Val exportReferencesData_ = AS3_Function( NULL, exportReferencesData );
	AS3_Val pushDataToRefObject_ = AS3_Function( NULL, pushDataToRefObject );
	AS3_Val setupImagePyramid_ = AS3_Function( NULL, setupImagePyramid );


	AS3_Val result = AS3_Object("setupGlobalBuffers: AS3ValType, setupImageHolders: AS3ValType, createRefObject: AS3ValType, pushImageToRefObject: AS3ValType, clearReferenceObjects: AS3ValType, buildRefIndex: AS3ValType, runTask: AS3ValType, runTask2: AS3ValType, getDataPointers: AS3ValType, clearBuffers: AS3ValType, exportReferencesData: AS3ValType, pushDataToRefObject: AS3ValType, setupImagePyramid: AS3ValType",
	setupGlobalBuffers_, setupImageHolders_, createRefObject_, pushImageToRefObject_, clearReferenceObjects_, buildRefIndex_, runTask_, runTask2_, getDataPointers_, clearBuffers_, exportReferencesData_, pushDataToRefObject_, setupImagePyramid_);

	AS3_Release( setupGlobalBuffers_ );
	AS3_Release( setupImageHolders_ );
	AS3_Release( createRefObject_ );
	AS3_Release( pushImageToRefObject_ );
	AS3_Release( clearReferenceObjects_ );
	AS3_Release( buildRefIndex_ );
	AS3_Release( runTask_ );
	AS3_Release( runTask2_ );
	AS3_Release( getDataPointers_ );
	AS3_Release( clearBuffers_ );
	AS3_Release( exportReferencesData_ );
	AS3_Release( pushDataToRefObject_ );
	AS3_Release( setupImagePyramid_ );

	AS3_LibInit( result );
	
	return 0;
}
