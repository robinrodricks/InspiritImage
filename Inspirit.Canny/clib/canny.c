#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "AS3.h"

static void dispose();
static int performDetection(register short *dx, register short *dy, register int *magn,
                       const long width, const long height, const int thresh,
                       const int thresh2, int *output);
static void ByteGradientFull(int *img, short *dx, short *dy, int *mag, int *max_mag);

int *image;
int *edges;

short *dx;
short *dy;
int *mag;

int **stktop;
int **stack;

int height = 240;
int width = 320;
int area = 0;

double thresh_low = 0.20;
double thresh_high = 0.60;

static int const EDGE_PIXEL = 16777215;
static int const NO_EDGE_PIXEL = 0;
static int const DUMMY_PIXEL = 2;

static long offset1[4];
static long offset2[4];
static int hist[200000];

static AS3_Val setupCanny(void* self, AS3_Val args)
{
	AS3_ArrayValue(args, "IntType, IntType, DoubleType, DoubleType", &width, &height, &thresh_low, &thresh_high);

	area = width * height;
	
	dispose();

	image = (int *) calloc(area, sizeof(int));
	edges = (int *) calloc(area, sizeof(int));

	dx = (short *) calloc(area, sizeof(short));
	dy = (short *) calloc(area, sizeof(short));
	mag = (int *) calloc(area, sizeof(int));

	stktop = stack = (int **)calloc(area, sizeof(int *));

	offset1[0] = width;	offset2[0] = width-1;    /* 225..270 and 45..90 */
	offset1[1] = width-1;	offset2[1] = -1;       /* 180..225 and 0..45 */
	offset1[2] = 1;		offset2[2] = width+1;    /* 315..360 and 135..180 */
	offset1[3] = width+1;	offset2[3] = width;      /* 270..315 and 90..135 */

	return 0;
}

static AS3_Val destroyCanny(void* self, AS3_Val args)
{
	dispose();
	return 0;
}

static AS3_Val getPointers(void* self, AS3_Val args)
{
	AS3_Val pointers = AS3_Array("AS3ValType", NULL);
	
	AS3_Set(pointers, AS3_Int(0), AS3_Ptr(&image));
	AS3_Set(pointers, AS3_Int(1), AS3_Ptr(&edges));
	AS3_Set(pointers, AS3_Int(2), AS3_Ptr(&thresh_low));
	AS3_Set(pointers, AS3_Int(3), AS3_Ptr(&thresh_high));
	
	return pointers;
}

static void dispose()
{
	if(image) free( image );
	if(edges) free( edges );
	if(dx) free( dx );
	if(dy) free( dy );
	if(mag) free( mag );
	if(stack) free( stack );
}

static AS3_Val runCanny(void* self, AS3_Val args)
{
	int i, numedges, highcount, maximum_mag, max_mag = 0, highthreshold, lowthreshold;
	register int *hp = hist;

	ByteGradientFull(image, dx, dy, mag, &max_mag);
	max_mag++;

	for(i = 0; i < max_mag; i++) *(hp++) = 0;

	register int *mg = mag;
	for(i = 0; i < area; i++)
	{
		hist[ (*(mg++)) ]++;
	}

	for(i = 1, numedges = 0, hp = hist+1; i < max_mag; i++, hp++)
	{
		if(*(hp) != 0) {
			maximum_mag = i;
			numedges += *(hp);
		}
	}

	highcount = (int)(numedges * thresh_high + 0.5);

	i = 1;
	numedges = hist[1];
	hp = hist+2;
	while( i<(maximum_mag-1) && numedges < highcount )
	{
		i++;
		numedges += *(hp++);
	}

	highthreshold = i;
	lowthreshold = (int)(highthreshold * thresh_low + 0.5);

	performDetection(dx, dy, mag, width, height, lowthreshold, highthreshold, edges);


	return 0;
}

inline int NMS(const int *magn, int mag, short dx, short dy, long w)
{
    int reg;
    long o1, o2;
    int i1,i2,sx,sy,s, m1,m2,denom;
    int answer;

    sx = dx < 0?-1:1;
    sy = dy < 0?-1:1;

    dx *= sx;
    dy *= sy;
    s = sx*sy;

    if (dy == 0) {
        m1 = *(magn + 1);
        m2 = *(magn - 1);
    } else if (dx == 0) {
        m1 = *(magn + w);
        m2 = *(magn - w);
    } else {
        if (s > 0) {
            if (dy <= dx) {
                reg = 2;
                i1 = dy;
                i2 = dx - dy;
                denom = dx;
            } else {
                reg = 3;
                i2 = dx;
                i1 = dy - dx;
                denom = dy;
            }
        } else { /* s < 0 */
            if (dy <= dx) {
                reg = 1;
                i2 = dy;
                i1 = dx - dy;
                denom = dx;
            } else {
                reg = 0;
                i1 = dx;
                i2 = dy - dx;
                denom = dy;
            }
        }
        o1 = offset1[reg];
        o2 = offset2[reg];
        m1 = (*(magn + o1)*i2 + *(magn + o2)*i1);
        m2 = (*(magn - o1)*i2 + *(magn - o2)*i1);
        mag *= denom;
    }
    answer = (mag>=m1 && mag>=m2 && m1!=m2);
    return answer; /* return 1 if passes NMS */
}

static int performDetection(register short *dx, register short *dy, register int *magn,
                       const long x_m, const long y_m, const int thresh,
                       const int thresh2, int *output)
{
    register int *out;
    int *magn_max, *magn_max_x;
    int mag, i;

    stack = stktop = 0;

    out  = output;
    while (out < output + x_m)               /* make first row 0 */
        *out++ = NO_EDGE_PIXEL;

    dx  += x_m + 1;                          /* skip over first col */
    dy  += x_m + 1;

    for (magn_max = magn + x_m*(y_m - 1) + 1, magn += x_m + 1; magn < magn_max; ) {
        *out++ = NO_EDGE_PIXEL;
        for (magn_max_x = magn + x_m - 2; magn < magn_max_x;
             out++, dx++, dy++, magn++) {
            if ((mag = *magn) < thresh) {            /* no possible edge */
                *out = NO_EDGE_PIXEL;
            } else {
                if (NMS(magn, mag, *dx, *dy, x_m)) { /* check if passes NMS */
                    if (mag >= thresh2) {
                        *out = EDGE_PIXEL;        /* definitely have an edge */
                        *(stktop++) = out;        /* put edge pixel addr on stack */
                    } else {
                        *out = DUMMY_PIXEL;       /* maybe have an edge */
                    }
                } else {
                    *out = NO_EDGE_PIXEL;         /* no edge here */
                }
            }
        }

        *out++ = NO_EDGE_PIXEL;

        dx   += 2;                  /* skip over last pixel on this line */
        dy   += 2;                  /* and first pixel on next line */
        magn += 2;
    }

    for (i=0; i<x_m; i++)
        *out++ = NO_EDGE_PIXEL;     /* Last row gets 0 too */

    while (stktop > stack) {
        out = *--stktop;

        /* look at neighbors. if a neighbor is DUMMY make it EDGE and add to stack */
        /* continue until stack empty */

        if (*(out -= x_m + 1) == DUMMY_PIXEL)			/* upper-left */
            *out = EDGE_PIXEL, *(stktop++) = out;
        if (*(++out) == DUMMY_PIXEL)				/* upper-middle */
            *out = EDGE_PIXEL, *(stktop++) = out;
        if (*(++out) == DUMMY_PIXEL)				/* upper-right */
            *out = EDGE_PIXEL, *(stktop++) = out;
        if (*(out += x_m) == DUMMY_PIXEL)				/* middle-right */
            *out = EDGE_PIXEL, *(stktop++) = out;
        if (*(out -= 2) == DUMMY_PIXEL)				/* middle-left */
            *out = EDGE_PIXEL, *(stktop++) = out;
        if (*(out += x_m) == DUMMY_PIXEL)				/* lower-left */
            *out = EDGE_PIXEL, *(stktop++) = out;
        if (*(++out) == DUMMY_PIXEL)				/* lower-middle */
            *out = EDGE_PIXEL, *(stktop++) = out;
        if (*(++out) == DUMMY_PIXEL)				/* lower-right */
            *out = EDGE_PIXEL, *(stktop++) = out;
    }

    {
        register int *outend = output + area; /* get rid of any remaining DUMMY's */

        for (out = output; out < outend; out++ )
            if (*out == DUMMY_PIXEL)
                *out = NO_EDGE_PIXEL;
    }

    return 0;
}

static void ByteGradientFull(int *img, short *dx, short *dy, int *mag, int *max_mag)
{
    register int *imgp, *img_xendp;
    register short *dxp, *dyp;
    register int *magp;
    register short x, y, a, b, c, d;
    register int scanpad;
    register int *img_endp, mmag;

    scanpad = width;

    imgp = img;
    dxp  = dx;
    dyp  = dy;
    magp = mag;

    for (img_endp = imgp + scanpad*(height-1); imgp < img_endp; ) {
        a = *imgp;
        c = *(imgp+scanpad);

        for (img_xendp = imgp + width - 1; imgp < img_xendp; ) {
            imgp++;

            b = *(imgp);
            d = *(imgp+scanpad);

            a = d - a;
            c = b - c;
            x = a + c;
            y = a - c;

            a = b;
            c = d;

            *(dxp++)  = x;
            *(dyp++)  = y;
		*(magp++) = (mmag = x*x + y*y);
            if(mmag > *max_mag) *max_mag = mmag;
        }

        imgp ++;
        *(dxp++) = *(dyp++) = *(magp++) = 0;
    }

    for (img_endp += width; imgp < img_endp; imgp++)
        *(dxp++) = *(dyp++) = *(magp++) = 0;
}

int main()
{
	AS3_Val runCanny_m = AS3_Function( NULL, runCanny );
	AS3_Val setupCanny_m = AS3_Function( NULL, setupCanny );
	AS3_Val getPointers_m = AS3_Function( NULL, getPointers );
	AS3_Val destroyCanny_m = AS3_Function( NULL, destroyCanny );

	AS3_Val result = AS3_Object(
						"runCanny: AS3ValType, setupCanny: AS3ValType, getPointers: AS3ValType, destroyCanny: AS3ValType",
						runCanny_m, setupCanny_m, getPointers_m, destroyCanny_m);

	AS3_Release( runCanny_m );
	AS3_Release( setupCanny_m );
	AS3_Release( getPointers_m );
	AS3_Release( destroyCanny_m );

	AS3_LibInit( result );

	return 0;
}
