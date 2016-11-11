#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>

#include "lkpyramid.h"

#define CV_8TO32F(x)  icv8x32fTab_cv[(x)+256]

#define  CV_MALLOC_ALIGN    32

inline void* icvAlignPtr( void* ptr, int align )
{
    return (void*)( ((size_t)ptr + align - 1) & -align );
}

inline size_t cvAlign( size_t size, int align ) 
{ 
    assert( (align & (align-1)) == 0 );
    return (size + align - 1) & ~(size_t)(align-1); 
}

void* cvAlloc( size_t size)
{
    uchar* udata = (uchar*)malloc(size + sizeof(void*) + CV_MALLOC_ALIGN);
    if(!udata) return 0;
    uchar** adata = icvAlignPtr((uchar**)udata + 1, CV_MALLOC_ALIGN);
    adata[-1] = udata;
    return adata;
}

void cvFree( void* ptr )
{
    if(ptr)
    {
        uchar* udata = ((uchar**)ptr)[-1];
        free(udata);
    }
}

static void
intersect( CvPoint2D32f pt, CvSize win_size, CvSize imgSize,
           CvPoint* min_pt, CvPoint* max_pt )
{
    CvPoint ipt;

    ipt.x = floor( pt.x );
    ipt.y = floor( pt.y );

    ipt.x -= win_size.width;
    ipt.y -= win_size.height;

    win_size.width = win_size.width * 2 + 1;
    win_size.height = win_size.height * 2 + 1;

    min_pt->x = IMAX( 0, -ipt.x );
    min_pt->y = IMAX( 0, -ipt.y );
    max_pt->x = IMIN( win_size.width, imgSize.width - ipt.x );
    max_pt->y = IMIN( win_size.height, imgSize.height - ipt.y );
}


/* compute dI/dx and dI/dy */
static void
icvCalcIxIy_32f( const double* src, int src_step, double* dstX, double* dstY, int dst_step,
                 CvSize src_size, const double* smooth_k, double* buffer0 )
{
    int src_width = src_size.width, dst_width = src_size.width-2;
    int x, height = src_size.height - 2;
    double* buffer1 = buffer0 + src_width;

    for( ; height--; src += src_step, dstX += dst_step, dstY += dst_step )
    {
        const double* src2 = src + src_step;
        const double* src3 = src + src_step*2;

        for( x = 0; x < src_width; x++ )
        {
            double t0 = (src3[x] + src[x])*smooth_k[0] + src2[x]*smooth_k[1];
            double t1 = src3[x] - src[x];
            buffer0[x] = t0; buffer1[x] = t1;
        }

        for( x = 0; x < dst_width; x++ )
        {
            double t0 = buffer0[x+2] - buffer0[x];
            double t1 = (buffer1[x] + buffer1[x+2])*smooth_k[0] + buffer1[x+1]*smooth_k[1];
            dstX[x] = t0; dstY[x] = t1;
        }
    }
}

static const double icv8x32fTab_cv[768] =
{
    -256.0, -255.0, -254.0, -253.0, -252.0, -251.0, -250.0, -249.0,
    -248.0, -247.0, -246.0, -245.0, -244.0, -243.0, -242.0, -241.0,
    -240.0, -239.0, -238.0, -237.0, -236.0, -235.0, -234.0, -233.0,
    -232.0, -231.0, -230.0, -229.0, -228.0, -227.0, -226.0, -225.0,
    -224.0, -223.0, -222.0, -221.0, -220.0, -219.0, -218.0, -217.0,
    -216.0, -215.0, -214.0, -213.0, -212.0, -211.0, -210.0, -209.0,
    -208.0, -207.0, -206.0, -205.0, -204.0, -203.0, -202.0, -201.0,
    -200.0, -199.0, -198.0, -197.0, -196.0, -195.0, -194.0, -193.0,
    -192.0, -191.0, -190.0, -189.0, -188.0, -187.0, -186.0, -185.0,
    -184.0, -183.0, -182.0, -181.0, -180.0, -179.0, -178.0, -177.0,
    -176.0, -175.0, -174.0, -173.0, -172.0, -171.0, -170.0, -169.0,
    -168.0, -167.0, -166.0, -165.0, -164.0, -163.0, -162.0, -161.0,
    -160.0, -159.0, -158.0, -157.0, -156.0, -155.0, -154.0, -153.0,
    -152.0, -151.0, -150.0, -149.0, -148.0, -147.0, -146.0, -145.0,
    -144.0, -143.0, -142.0, -141.0, -140.0, -139.0, -138.0, -137.0,
    -136.0, -135.0, -134.0, -133.0, -132.0, -131.0, -130.0, -129.0,
    -128.0, -127.0, -126.0, -125.0, -124.0, -123.0, -122.0, -121.0,
    -120.0, -119.0, -118.0, -117.0, -116.0, -115.0, -114.0, -113.0,
    -112.0, -111.0, -110.0, -109.0, -108.0, -107.0, -106.0, -105.0,
    -104.0, -103.0, -102.0, -101.0, -100.0,  -99.0,  -98.0,  -97.0,
     -96.0,  -95.0,  -94.0,  -93.0,  -92.0,  -91.0,  -90.0,  -89.0,
     -88.0,  -87.0,  -86.0,  -85.0,  -84.0,  -83.0,  -82.0,  -81.0,
     -80.0,  -79.0,  -78.0,  -77.0,  -76.0,  -75.0,  -74.0,  -73.0,
     -72.0,  -71.0,  -70.0,  -69.0,  -68.0,  -67.0,  -66.0,  -65.0,
     -64.0,  -63.0,  -62.0,  -61.0,  -60.0,  -59.0,  -58.0,  -57.0,
     -56.0,  -55.0,  -54.0,  -53.0,  -52.0,  -51.0,  -50.0,  -49.0,
     -48.0,  -47.0,  -46.0,  -45.0,  -44.0,  -43.0,  -42.0,  -41.0,
     -40.0,  -39.0,  -38.0,  -37.0,  -36.0,  -35.0,  -34.0,  -33.0,
     -32.0,  -31.0,  -30.0,  -29.0,  -28.0,  -27.0,  -26.0,  -25.0,
     -24.0,  -23.0,  -22.0,  -21.0,  -20.0,  -19.0,  -18.0,  -17.0,
     -16.0,  -15.0,  -14.0,  -13.0,  -12.0,  -11.0,  -10.0,   -9.0,
      -8.0,   -7.0,   -6.0,   -5.0,   -4.0,   -3.0,   -2.0,   -1.0,
       0.0,    1.0,    2.0,    3.0,    4.0,    5.0,    6.0,    7.0,
       8.0,    9.0,   10.0,   11.0,   12.0,   13.0,   14.0,   15.0,
      16.0,   17.0,   18.0,   19.0,   20.0,   21.0,   22.0,   23.0,
      24.0,   25.0,   26.0,   27.0,   28.0,   29.0,   30.0,   31.0,
      32.0,   33.0,   34.0,   35.0,   36.0,   37.0,   38.0,   39.0,
      40.0,   41.0,   42.0,   43.0,   44.0,   45.0,   46.0,   47.0,
      48.0,   49.0,   50.0,   51.0,   52.0,   53.0,   54.0,   55.0,
      56.0,   57.0,   58.0,   59.0,   60.0,   61.0,   62.0,   63.0,
      64.0,   65.0,   66.0,   67.0,   68.0,   69.0,   70.0,   71.0,
      72.0,   73.0,   74.0,   75.0,   76.0,   77.0,   78.0,   79.0,
      80.0,   81.0,   82.0,   83.0,   84.0,   85.0,   86.0,   87.0,
      88.0,   89.0,   90.0,   91.0,   92.0,   93.0,   94.0,   95.0,
      96.0,   97.0,   98.0,   99.0,  100.0,  101.0,  102.0,  103.0,
     104.0,  105.0,  106.0,  107.0,  108.0,  109.0,  110.0,  111.0,
     112.0,  113.0,  114.0,  115.0,  116.0,  117.0,  118.0,  119.0,
     120.0,  121.0,  122.0,  123.0,  124.0,  125.0,  126.0,  127.0,
     128.0,  129.0,  130.0,  131.0,  132.0,  133.0,  134.0,  135.0,
     136.0,  137.0,  138.0,  139.0,  140.0,  141.0,  142.0,  143.0,
     144.0,  145.0,  146.0,  147.0,  148.0,  149.0,  150.0,  151.0,
     152.0,  153.0,  154.0,  155.0,  156.0,  157.0,  158.0,  159.0,
     160.0,  161.0,  162.0,  163.0,  164.0,  165.0,  166.0,  167.0,
     168.0,  169.0,  170.0,  171.0,  172.0,  173.0,  174.0,  175.0,
     176.0,  177.0,  178.0,  179.0,  180.0,  181.0,  182.0,  183.0,
     184.0,  185.0,  186.0,  187.0,  188.0,  189.0,  190.0,  191.0,
     192.0,  193.0,  194.0,  195.0,  196.0,  197.0,  198.0,  199.0,
     200.0,  201.0,  202.0,  203.0,  204.0,  205.0,  206.0,  207.0,
     208.0,  209.0,  210.0,  211.0,  212.0,  213.0,  214.0,  215.0,
     216.0,  217.0,  218.0,  219.0,  220.0,  221.0,  222.0,  223.0,
     224.0,  225.0,  226.0,  227.0,  228.0,  229.0,  230.0,  231.0,
     232.0,  233.0,  234.0,  235.0,  236.0,  237.0,  238.0,  239.0,
     240.0,  241.0,  242.0,  243.0,  244.0,  245.0,  246.0,  247.0,
     248.0,  249.0,  250.0,  251.0,  252.0,  253.0,  254.0,  255.0,
     256.0,  257.0,  258.0,  259.0,  260.0,  261.0,  262.0,  263.0,
     264.0,  265.0,  266.0,  267.0,  268.0,  269.0,  270.0,  271.0,
     272.0,  273.0,  274.0,  275.0,  276.0,  277.0,  278.0,  279.0,
     280.0,  281.0,  282.0,  283.0,  284.0,  285.0,  286.0,  287.0,
     288.0,  289.0,  290.0,  291.0,  292.0,  293.0,  294.0,  295.0,
     296.0,  297.0,  298.0,  299.0,  300.0,  301.0,  302.0,  303.0,
     304.0,  305.0,  306.0,  307.0,  308.0,  309.0,  310.0,  311.0,
     312.0,  313.0,  314.0,  315.0,  316.0,  317.0,  318.0,  319.0,
     320.0,  321.0,  322.0,  323.0,  324.0,  325.0,  326.0,  327.0,
     328.0,  329.0,  330.0,  331.0,  332.0,  333.0,  334.0,  335.0,
     336.0,  337.0,  338.0,  339.0,  340.0,  341.0,  342.0,  343.0,
     344.0,  345.0,  346.0,  347.0,  348.0,  349.0,  350.0,  351.0,
     352.0,  353.0,  354.0,  355.0,  356.0,  357.0,  358.0,  359.0,
     360.0,  361.0,  362.0,  363.0,  364.0,  365.0,  366.0,  367.0,
     368.0,  369.0,  370.0,  371.0,  372.0,  373.0,  374.0,  375.0,
     376.0,  377.0,  378.0,  379.0,  380.0,  381.0,  382.0,  383.0,
     384.0,  385.0,  386.0,  387.0,  388.0,  389.0,  390.0,  391.0,
     392.0,  393.0,  394.0,  395.0,  396.0,  397.0,  398.0,  399.0,
     400.0,  401.0,  402.0,  403.0,  404.0,  405.0,  406.0,  407.0,
     408.0,  409.0,  410.0,  411.0,  412.0,  413.0,  414.0,  415.0,
     416.0,  417.0,  418.0,  419.0,  420.0,  421.0,  422.0,  423.0,
     424.0,  425.0,  426.0,  427.0,  428.0,  429.0,  430.0,  431.0,
     432.0,  433.0,  434.0,  435.0,  436.0,  437.0,  438.0,  439.0,
     440.0,  441.0,  442.0,  443.0,  444.0,  445.0,  446.0,  447.0,
     448.0,  449.0,  450.0,  451.0,  452.0,  453.0,  454.0,  455.0,
     456.0,  457.0,  458.0,  459.0,  460.0,  461.0,  462.0,  463.0,
     464.0,  465.0,  466.0,  467.0,  468.0,  469.0,  470.0,  471.0,
     472.0,  473.0,  474.0,  475.0,  476.0,  477.0,  478.0,  479.0,
     480.0,  481.0,  482.0,  483.0,  484.0,  485.0,  486.0,  487.0,
     488.0,  489.0,  490.0,  491.0,  492.0,  493.0,  494.0,  495.0,
     496.0,  497.0,  498.0,  499.0,  500.0,  501.0,  502.0,  503.0,
     504.0,  505.0,  506.0,  507.0,  508.0,  509.0,  510.0,  511.0
};




static const void*
icvAdjustRect( const void* srcptr, int src_step, int pix_size,
               CvSize src_size, CvSize win_size,
               CvPoint ip, CvRect* pRect )
{
    CvRect rect;
    const char* src = (const char*)srcptr;

    if( ip.x >= 0 )
    {
        src += ip.x*pix_size;
        rect.x = 0;
    }
    else
    {
        rect.x = -ip.x;
        if( rect.x > win_size.width )
            rect.x = win_size.width;
    }

    if( ip.x + win_size.width < src_size.width )
        rect.width = win_size.width;
    else
    {
        rect.width = src_size.width - ip.x - 1;
        if( rect.width < 0 )
        {
            src += rect.width*pix_size;
            rect.width = 0;
        }
        assert( rect.width <= win_size.width );
    }

    if( ip.y >= 0 )
    {
        src += ip.y * src_step;
        rect.y = 0;
    }
    else
        rect.y = -ip.y;

    if( ip.y + win_size.height < src_size.height )
        rect.height = win_size.height;
    else
    {
        rect.height = src_size.height - ip.y - 1;
        if( rect.height < 0 )
        {
            src += rect.height*src_step;
            rect.height = 0;
        }
    }

	pRect->x = rect.x;
	pRect->y = rect.y;
    pRect->width = rect.width;
	pRect->height = rect.height;
	
    return src - rect.x*pix_size;
}

static int  mycvGetRectSubPix_8u32f_C1R
( const uchar* src, int src_step, CvSize src_size,
  double* dst, int dst_step, CvSize win_size, CvPoint2D32f center )
{
    CvPoint ip;
    double  a12, a22, b1, b2;
    double a, b;
    double s = 0;
    int i, j;

    center.x -= (double)(win_size.width-1)*(double)0.5;
    center.y -= (double)(win_size.height-1)*(double)0.5;

    ip.x = floor( center.x );
    ip.y = floor( center.y );

    if( win_size.width <= 0 || win_size.height <= 0 ) return -1;//CV_BADRANGE_ERR;

    a = center.x - (double)ip.x;
    b = center.y - (double)ip.y;
    a = FMAX(a,0.0001);
    a12 = a*((double)1.0-b);
    a22 = a*b;
    b1 = (double)1.0 - b;
    b2 = b;
    s = ((double)1.0 - a)/a;

    if( 0 <= ip.x && ip.x + win_size.width < src_size.width &&
        0 <= ip.y && ip.y + win_size.height < src_size.height )
    {
        // extracted rectangle is totally inside the image
        src += ip.y * src_step + ip.x;

        for( ; win_size.height--; src += src_step, dst += dst_step )
        {
            double prev = ((double)1. - a)*(b1*CV_8TO32F(src[0]) + b2*CV_8TO32F(src[src_step]));
            for( j = 0; j < win_size.width; j++ )
            {
                double t = a12*CV_8TO32F(src[j+1]) + a22*CV_8TO32F(src[j+1+src_step]);
                dst[j] = prev + t;
                prev = (double)(t*s);
            }
        }
    }
    else
    {
        CvRect r;

        src = (const uchar*)icvAdjustRect( src, src_step, 1, src_size, win_size, ip, &r);

        for( i = 0; i < win_size.height; i++, dst += dst_step )
        {
            const uchar *src2 = src + src_step;

            if( i < r.y || i >= r.height )
                src2 -= src_step;

            for( j = 0; j < r.x; j++ )
            {
                double s0 = CV_8TO32F(src[r.x])*b1 +
                           CV_8TO32F(src2[r.x])*b2;

                dst[j] = (double)(s0);
            }

            if( j < r.width )
            {
                double prev = ((double)1. - a)*(b1*CV_8TO32F(src[j]) + b2*CV_8TO32F(src2[j]));

                for( ; j < r.width; j++ )
                {
                    double t = a12*CV_8TO32F(src[j+1]) + a22*CV_8TO32F(src2[j+1]);
                    dst[j] = prev + t;
                    prev = (double)(t*s);
                }
            }

            for( ; j < win_size.width; j++ )
            {
                double s0 = CV_8TO32F(src[r.width])*b1 +
                           CV_8TO32F(src2[r.width])*b2;

                dst[j] = (double)(s0);
            }

            if( i < r.height )
                src = src2;
        }
    }

    return 0;
}

uchar *pyrPatchBuffer = 0;

int    *step = 0;
double *scale = 0;
CvSize *size = 0;

double* patchI;
double* patchJ;
double* Ix;
double* Iy;

void initOpticalFlow( uchar* arrA[3], uchar* arrB[3], 
						const int sizeW[3], const int sizeH[3], int patch_size, int flags)
{
	free( step );
	free( scale );
	free( size );
    free( pyrPatchBuffer );
	
	int level = 3;
    int bufferBytes = 0;
    int i;

    step = (int *)malloc( level * sizeof(int) );
    scale = (double *)malloc( level * sizeof(double) );
    size = (CvSize *)malloc( level * sizeof(CvSize) );

    for (i = 0; i < level; i++) 
	{
	    step[i] = sizeW[i];
		CvSize *lsz = &size[i];
		lsz->width = sizeW[i];
		lsz->height = sizeH[i];
	    scale[i] = (double)1.0 / (double)(1<<i);
    }
	
	int patchSize = patch_size * 2 + 1;
    int patchLen = patchSize * patchSize;
    int srcPatchLen = (patchSize + 2)*(patchSize + 2);
	
	bufferBytes = (srcPatchLen + patchLen * 3) * sizeof( double );
	pyrPatchBuffer = (uchar*)malloc( bufferBytes );

	patchI = (double*)pyrPatchBuffer;
	patchJ = patchI + srcPatchLen;
	Ix = patchJ + patchLen;
	Iy = Ix + patchLen;
}



void myCalcOpticalFlowPyrLK( uchar* arrA[3], uchar* arrB[3],
						const int sizeW[3], const int sizeH[3],
                        const CvPoint2D32f * featuresA, CvPoint2D32f * featuresB,
                        int count, int patch_size, char *status, double *error,
                        int flags )
{    
	
	CvSize winSize;
	winSize.width = winSize.height = patch_size;
	
	CvTermCriteria criteria;
	criteria.type = CV_TERMCRIT_ITER|CV_TERMCRIT_EPS;
	criteria.max_iter = 3;
	criteria.epsilon = 0.1 * 0.1;

    const double smoothKernel[3] = { 0.09375, 0.3125, 0.09375 };  /* 3/32, 10/32, 3/32 */   

    int i, l;

    CvSize patchSize;
	patchSize.width = winSize.width * 2 + 1;
	patchSize.height = winSize.height * 2 + 1;

	CvSize imgSize;
    imgSize.width = sizeW[0];
    imgSize.height = sizeH[0];

    if( count == 0 ) return;
    
    int level = 3 - 1;

	char *_status;
    if( !status ) ( status = _status = (char*)malloc( count*sizeof(char) ));    

    memset( status, 1, count );
    if( error ) memset( error, 0, count*sizeof(double) );

    if( !(flags & CV_LKFLOW_INITIAL_GUESSES) ) memcpy( featuresB, featuresA, count*sizeof(CvPoint2D32f));

    /* do processing from top pyramid level (smallest image)
       to the bottom (original image) */
    for( l = level; l >= 0; l-- )
    {
        CvSize levelSize = size[l];
        int levelStep = step[l];

        {

        /* find flow for each given point */
        for( i = 0; i < count; i++ )
        {
            CvPoint2D32f v;
            CvPoint minI, maxI, minJ, maxJ;
            CvSize isz, jsz;
            int pt_status;
            CvPoint2D32f u;
            CvPoint prev_minJ, prev_maxJ;
			prev_minJ.x = -1;
			prev_minJ.y = -1;
			prev_maxJ.x = -1;
			prev_maxJ.y = -1;
            double Gxx = 0, Gxy = 0, Gyy = 0, D = 0, minEig = 0;
            double prev_mx = 0, prev_my = 0;
            int j, x, y;

            v.x = featuresB[i].x;
            v.y = featuresB[i].y;
            if( l < level )
            {
                v.x += v.x;
                v.y += v.y;
            }
            else
            {
                v.x = (double)(v.x * scale[l]);
                v.y = (double)(v.y * scale[l]);
            }

            pt_status = status[i];
            if( !pt_status ) continue;

            minI.x = maxI.x = minJ.x = maxJ.x = 0;
            minI.y = maxI.y = minJ.y = maxJ.y = 0;

            u.x = (double) (featuresA[i].x * scale[l]);
            u.y = (double) (featuresA[i].y * scale[l]);

            intersect( u, winSize, levelSize, &minI, &maxI );
            isz.width = jsz.width = maxI.x - minI.x + 2;
			isz.height = jsz.height = maxI.y - minI.y + 2;
            u.x += (minI.x - (patchSize.width - maxI.x + 1))*(double)0.5;
            u.y += (minI.y - (patchSize.height - maxI.y + 1))*(double)0.5;

            if( isz.width < 3 || isz.height < 3 ||
                mycvGetRectSubPix_8u32f_C1R( arrA[l], levelStep, levelSize, patchI, isz.width, isz, u ) < 0 )
            {
                /* point is outside the image. take the next */
                status[i] = 0;
                continue;
            }

            icvCalcIxIy_32f( patchI, isz.width, Ix, Iy, (isz.width-2), isz, smoothKernel, patchJ );

            for( j = 0; j < criteria.max_iter; j++ )
            {
                double bx = 0, by = 0;
                double mx, my;
                CvPoint2D32f _v;

                intersect( v, winSize, levelSize, &minJ, &maxJ );

                minJ.x = IMAX( minJ.x, minI.x );
                minJ.y = IMAX( minJ.y, minI.y );

                maxJ.x = IMIN( maxJ.x, maxI.x );
                maxJ.y = IMIN( maxJ.y, maxI.y );

                jsz.width = maxJ.x - minJ.x;
				jsz.height = maxJ.y - minJ.y;

                _v.x = v.x + (minJ.x - (patchSize.width - maxJ.x + 1))*(double)0.5;
                _v.y = v.y + (minJ.y - (patchSize.height - maxJ.y + 1))*(double)0.5;

                if( jsz.width < 1 || jsz.height < 1 ||
                    mycvGetRectSubPix_8u32f_C1R( arrB[l], levelStep, levelSize, patchJ, jsz.width, jsz, _v ) < 0 )
                {
                    /* point is outside image. take the next */
                    pt_status = 0;
                    break;
                }

                if( maxJ.x == prev_maxJ.x && maxJ.y == prev_maxJ.y &&
                    minJ.x == prev_minJ.x && minJ.y == prev_minJ.y )
                {
                    for( y = 0; y < jsz.height; y++ )
                    {
                        const double* pi = patchI +
                            (y + minJ.y - minI.y + 1)*isz.width + minJ.x - minI.x + 1;
                        const double* pj = patchJ + y*jsz.width;
                        const double* ix = Ix +
                            (y + minJ.y - minI.y)*(isz.width-2) + minJ.x - minI.x;
                        const double* iy = Iy + (ix - Ix);

                        for( x = 0; x < jsz.width; x++ )
                        {
                            double t0 = pi[x] - pj[x];
                            bx += t0 * ix[x];
                            by += t0 * iy[x];
                        }
                    }
                }
                else
                {
                    Gxx = Gyy = Gxy = 0;
                    for( y = 0; y < jsz.height; y++ )
                    {
                        const double* pi = patchI +
                            (y + minJ.y - minI.y + 1)*isz.width + minJ.x - minI.x + 1;
                        const double* pj = patchJ + y*jsz.width;
                        const double* ix = Ix +
                            (y + minJ.y - minI.y)*(isz.width-2) + minJ.x - minI.x;
                        const double* iy = Iy + (ix - Ix);

                        for( x = 0; x < jsz.width; x++ )
                        {
                            double t = pi[x] - pj[x];
                            bx += (double) (t * ix[x]);
                            by += (double) (t * iy[x]);
                            Gxx += ix[x] * ix[x];
                            Gxy += ix[x] * iy[x];
                            Gyy += iy[x] * iy[x];
                        }
                    }

                    D = Gxx * Gyy - Gxy * Gxy;
                    if( D < DBL_EPSILON )
                    {
                        pt_status = 0;
                        break;
                    }

                    // Adi Shavit - 2008.05
                    if( flags & CV_LKFLOW_GET_MIN_EIGENVALS )
					{
                        minEig = (Gyy + Gxx - fast_sqrt((Gxx-Gyy)*(Gxx-Gyy) + 4.*Gxy*Gxy))/(2*jsz.height*jsz.width);
					}

                    D = (double)1. / D;

                    prev_minJ.x = minJ.x;
					prev_minJ.y = minJ.y;
                    prev_maxJ.x = maxJ.x;
					prev_maxJ.y = maxJ.y;
                }

                mx = (double) ((Gyy * bx - Gxy * by) * D);
                my = (double) ((Gxx * by - Gxy * bx) * D);

                v.x += mx;
                v.y += my;

                if( mx * mx + my * my < criteria.epsilon )
                    break;

                if( j > 0 && fabs(mx + prev_mx) < 0.01 && fabs(my + prev_my) < 0.01 )
                {
                    v.x -= mx*(double)0.5;
                    v.y -= my*(double)0.5;
                    break;
                }
                prev_mx = mx;
                prev_my = my;
            }

            featuresB[i] = v;
            status[i] = (char)pt_status;
            if( l == 0 && error && pt_status )
            {
                /* calc error */
                double err = 0;
                if( flags & CV_LKFLOW_GET_MIN_EIGENVALS )
                    err = minEig;
                else
                {
                    for( y = 0; y < jsz.height; y++ )
                    {
                        const double* pi = patchI +
                            (y + minJ.y - minI.y + 1)*isz.width + minJ.x - minI.x + 1;
                        const double* pj = patchJ + y*jsz.width;

                        for( x = 0; x < jsz.width; x++ )
                        {
                            double t = pi[x] - pj[x];
                            err += t * t;
                        }
                    }
                    err = fast_sqrt(err);
                }
                error[i] = (double)err;
            }
        } // end of point processing loop (i)
        }
    } // end of pyramid levels loop (l)
	
	if(_status) free( _status );
}
