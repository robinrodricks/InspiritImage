#include <math.h>
#include <time.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "AS3.h"

void linearSolver(int b, double *x, double *x0, double a, double c);
void linearSolverRGB(double a, double c);
void linearSolverUV(double a, double c);
void destroy();
void reset();
void drawFluidImage();
void fadeR();
void fadeRGB();
void advect(int b, double *_d, double *d0, double *du, double *dv);
void advectRGB(double *du, double *dv);
void project(double *x, double *y, double *p, double *div);
void setBoundary(int bound, double *x);
void setBoundaryRGB();
void drawLine(int *layer, int x0, int y0, int x1, int y1, int c);
int compareParticles(const void *vp1, const void *vp2);


static const double FLUID_DEFAULT_DT					= 0.5;
static const double FLUID_DEFAULT_VISC					= 0.0001;
static const double FLUID_DEFAULT_COLOR_DIFFUSION 		= 0.0;
static const double FLUID_DEFAULT_FADESPEED				= 0.003;
static const int FLUID_DEFAULT_SOLVER_ITERATIONS		= 10;
static const int FLUID_DEFAULT_VORTICITY_CONFINEMENT 	= 0;

static const int PARTICLES_MAX							= 25000;
static const double MOMENTUM 							= 0.5;
static const double FLUID_FORCE 						= 0.6;

int NX;
int NY;
int NX2;
int NY2;
int numCells;
int screenW;
int screenH;
double isw;
double ish;
double invNX;
double invNY;
double invNumCells;

double _dt;
int _isRGB = 1;
int _solverIterations;
double _colorDiffusion;
int _doVorticityConfinement = 0;

int wrap_x = 0;
int wrap_y = 0;

double _visc;
double _fadeSpeed;

double _avgDensity;			// this will hold the average color of the last frame (how full it is)
double _uniformity;			// this will hold the uniformity of the last frame (how uniform the color is);
double _avgSpeed;

double *r;
double *g;
double *b;
double *rOld;
double *gOld;
double *bOld;
double *u;
double *uOld;
double *v;
double *vOld;
double *curl_abs;
double *curl_orig;

int *fluidsImage;
int *particlesImage;

double *particlesPool;

double *particles;
double *particles2;
double *_particles;
double *_particles2;
int particlesNum = 0;
int drawMode = 0;

inline double FMIN(register double a, register double b)
{
	return ((a) < (b) ? (a) : (b));
}
inline int IMIN(register int a, register int b)
{
	return ((a) < (b) ? (a) : (b));
}
inline double FMAX(register double a, register double b)
{
	return ((a) < (b) ? (b) : (a));
}
inline int IMAX(register int a, register int b)
{
	return ((a) < (b) ? (b) : (a));
}
inline int FLUID_IX(register int i, register int j)
{
	return (int)((i) + ((NX2) * (j)));
}

inline void addSourceUV()
{
	register double *up, *vp, *uop, *vop;

	for(up = u, vp = v, uop = uOld, vop = vOld; up < u+numCells;)
	{
		*(up++) += (_dt * (*(uop++)));
		*(vp++) += (_dt * (*(vop++)));
	}
}

inline void addSourceRGB()
{
	register double *rp, *gp, *bp, *rop, *bop, *gop;
	for(rp = r, gp = g, bp = b, rop = rOld, bop = bOld, gop = gOld; rp < r+numCells;)
	{
		*(rp++) += (_dt * (*(rop++)));
		*(gp++) += (_dt * (*(gop++)));
		*(bp++) += (_dt * (*(bop++)));
	}
}

inline void addSource(register double *x, register double *x0)
{
	for(; x < x+numCells;)
	{
		*(x++) += (_dt * (*(x0++)));
	}
}

inline void diffuse(const int b, double *c, double *c0, const double _diff)
{
	const double a = _dt * _diff * NX * NY;
	linearSolver(b, c, c0, a, 1.0 + 4 * a);
}

inline void diffuseRGB(const double _diff)
{
	const double a = _dt * _diff * NX * NY;
	linearSolverRGB(a, 1.0 + 4 * a);
}

inline void diffuseUV(const double _diff)
{
	const double a = _dt * _diff * NX * NY;
	linearSolverUV(a, 1.0 + 4 * a);
}

inline void swapUV()
{
	register double *_tmp;
	_tmp = u;
	u = uOld;
	uOld = _tmp;

	_tmp = v;
	v = vOld;
	vOld = _tmp;
}

inline void swapR()
{
	register double *_tmp;
	_tmp = r;
	r = rOld;
	rOld = _tmp;
}

inline void swapRGB()
{
	register double *_tmp;
	_tmp = r;
	r = rOld;
	rOld = _tmp;

	_tmp = g;
	g = gOld;
	gOld = _tmp;

	_tmp = b;
	b = bOld;
	bOld = _tmp;
}

inline void drawLine(int *layer, int x0, int y0, int x1, int y1, const int c)
{
	int dy = y1 - y0;
	int dx = x1 - x0;
	int stepx, stepy, fraction;

	if (dy < 0) { dy = -dy;  stepy = -screenW; } else { stepy = screenW; }
	if (dx < 0) { dx = -dx;  stepx = -1; } else { stepx = 1; }
	dy <<= 1;
	dx <<= 1;

	y0 *= screenW;
	y1 *= screenW;
	layer[x0+y0] = c;
	if (dx > dy) {
		fraction = dy - (dx>>1);
		while (x0 != x1) {
			if (fraction >= 0) {
				y0 += stepy;
				fraction -= dx;
			}
			x0 += stepx;
			fraction += dy;
			layer[x0+y0] = c;
		}
	} else {
		fraction = dx - (dy>>1);
		while (y0 != y1) {
			if (fraction >= 0) {
				x0 += stepx;
				fraction -= dy;
			}
			y0 += stepy;
			fraction += dx;
			layer[x0+y0] = c;
		}
	}
}

void destroy()
{
	if(fluidsImage)	free(fluidsImage);
	if(particles)	free(particles);
	if(particlesImage) free(particlesImage);
	if(particlesPool)	free(particlesPool);

	if(particles2)	free(particles2);

    if(r)		free(r);
    if(rOld)	free(rOld);

    if(g)		free(g);
    if(gOld)	free(gOld);

    if(b)		free(b);
    if(bOld)	free(bOld);

    if(u)		free(u);
    if(uOld)	free(uOld);
    if(v)		free(v);
    if(vOld)	free(vOld);
    if(curl_abs)	free(curl_abs);
    if(curl_orig)	free(curl_orig);
}


void reset()
{
	destroy();

	r    = (double*)calloc( numCells, sizeof(double) );
	rOld = (double*)calloc( numCells, sizeof(double) );

	g    = (double*)calloc( numCells, sizeof(double) );
	gOld = (double*)calloc( numCells, sizeof(double) );

	b    = (double*)calloc( numCells, sizeof(double) );
	bOld = (double*)calloc( numCells, sizeof(double) );

	u    = (double*)calloc( numCells, sizeof(double) );
	uOld = (double*)calloc( numCells, sizeof(double) );
	v    = (double*)calloc( numCells, sizeof(double) );
	vOld = (double*)calloc( numCells, sizeof(double) );

	curl_abs = (double*)calloc( numCells, sizeof(double) );
	curl_orig = (double*)calloc( numCells, sizeof(double) );

	fluidsImage = (int*)calloc( NX*NY, sizeof(int) );

	particlesImage = (int*)calloc( screenW*screenH, sizeof(int) );
	particles = (double*)calloc( PARTICLES_MAX*6, sizeof(double) );
	particlesPool = (double*)calloc( 50, sizeof(double) );

	particles2 = (double*)calloc( PARTICLES_MAX*6, sizeof(double) );

	_particles = particles;
	_particles2 = particles2;
	particlesNum = 0;

	srand( (unsigned)time(NULL) );
}

void drawFluidImage()
{
	int i, j;
	register int *ip;
	register double *rp, *gp, *bp;

	ip = fluidsImage;
	rp = r+NX2+1;
	gp = g+NX2+1;
	bp = b+NX2+1;

	for( i = 1; i < NY+1; i++ )
	{
		for( j = 1; j < NX+1; j++ )
		{
			*(ip++) = (int)( (int)((int)(*(rp++) * 0xFF)<<16) | (int)((int)(*(gp++) * 0xFF)<<8) | (int)(*(bp++) * 0xFF) );
		}
		rp+=2;
		gp+=2;
		bp+=2;
	}
}

void drawParticleImage()
{
	register double *pp;
	register double *pp2;
	int fluidIndex;
	double x, y, vx, vy, alpha, mass;

	int cnt = 0;

	for( pp = _particles, pp2 = _particles2; pp < _particles+particlesNum*6; )
	{
		alpha = *(pp++);
		x = *(pp++);
		y = *(pp++);
		vx = *(pp++);
		vy = *(pp++);
		mass = *(pp++);

		fluidIndex = (int)(x*isw*NX+1.5) + (int)(y*ish*NY+1.5) * NX2;
		vx = u[fluidIndex] * screenW * mass * FLUID_FORCE + vx * MOMENTUM;
		vy = v[fluidIndex] * screenH * mass * FLUID_FORCE + vy * MOMENTUM;

		x += vx;
		y += vy;

		if (x < 1.) {
			if (wrap_x == 1) {
				x = screenW-1;
			} else {
				x = 1.0;
				vx *= -1.0;
			}
		}
		else if (x > screenW) {
			if (wrap_x == 1) {
				x = 1;
			} else {
				x = screenW - 1;
				vx *= -1.0;
			}
		}

		if (y < 1.) {
			if (wrap_y == 1) {
				y = screenH - 1;
			} else {
				y = 1.0;
				vy *= -1.0;
			}
		}
		else if (y > screenH) {
			if (wrap_y == 1) {
				y = 1;
			} else {
				y = screenH - 1;
				vy *= -1.0;
			}
		}

		// hackish way to make particles glitter when they slow down a lot
		if(vx * vx + vy * vy < 1.0) {
			vx = (double)(((double)(rand()%201)-100.0) / 100.0);
			vy = (double)(((double)(rand()%201)-100.0) / 100.0);
		}

		/*a = (int)(alpha * 0xFF + 0.5);
		if(drawMode == 1 || drawMode == 2) {
			pix = (a<<24) | (a<<16) | (a<<8) | a;
		} else if(drawMode == 3 || drawMode == 4){
			pix = (0xFF<<24) | ((int)(r[fluidIndex]*0xFF)<<16) | ((int)(g[fluidIndex]*0xFF)<<8) | (int)(b[fluidIndex]*0xFF);
		}*/

		alpha = alpha * 0.996;

		if(alpha > 0.1){
			*(pp2++) = alpha;
			*(pp2++) = x;
			*(pp2++) = y;
			*(pp2++) = vx;
			*(pp2++) = vy;
			*(pp2++) = mass;
			cnt++;
		}

		//drawLine(particlesImage, (int)px, (int)py, (int)(x), (int)(y), pix);
	}

	particlesNum = cnt;
	register double *tmp = _particles;
	_particles = _particles2;
	_particles2 = tmp;
}

void addParticles(double x, double y, int num)
{
	register double *pp;

	if(particlesNum+num > PARTICLES_MAX) {
		for( pp = _particles; pp < _particles+(num)*6; )
		{
			*(pp++) = (double)((double)((rand()%(71))+30.0) / 100.0);
			*(pp++) = x + (double)((rand()%(21))-10);
			*(pp++) = y + (double)((rand()%(21))-10);
			*(pp++) = 0.0;
			*(pp++) = 0.0;
			*(pp++) = (double)((double)((rand()%(91))+10.0) / 100.0);
		}
	} else {
		for( pp = _particles+particlesNum*6; pp < _particles+(particlesNum+num)*6; )
		{
			*(pp++) = (double)((double)((rand()%(71))+30.0) / 100.0);
			*(pp++) = x + (double)((rand()%(21))-10);
			*(pp++) = y + (double)((rand()%(21))-10);
			*(pp++) = 0.0;
			*(pp++) = 0.0;
			*(pp++) = (double)((double)((rand()%(91))+10.0) / 100.0);
		}
		particlesNum += num;
	}
}

void calcVorticityConfinement(register double *_x, register double *_y)
{
	double dw_dx, dw_dy, length, vv;
	int i, j, index;

	// Calculate magnitude of (u,v) for each cell. (|w|)
	for (j = NY; j > 0; --j)
	{
		  index = FLUID_IX(NX, j);
		  for (i = NX; i > 0; --i)
		  {
			    dw_dy = u[index + NX2] - u[index - NX2];
			    dw_dx = v[index + 1] - v[index - 1];

			    vv = (dw_dy - dw_dx) * .5;

			    curl_orig[ index ] = vv;
			    curl_abs[ index ] = vv < 0 ? -vv : vv;

			    --index;
		  }
	}

	for (j = NY-1; j > 1; --j)
	{
		  index = FLUID_IX(NX-1, j);
		  for (i = NX-1; i > 1; --i)
		  {
			    dw_dx = curl_abs[index + 1] - curl_abs[index - 1];
			    dw_dy = curl_abs[index + NX2] - curl_abs[index - NX2];

			    length = sqrt(dw_dx * dw_dx + dw_dy * dw_dy) + 0.000001;

			    length = 2 / length;
			    dw_dx *= length;
			    dw_dy *= length;

			    vv = curl_orig[ index ];

			    _x[ index ] = dw_dy * -vv;
			    _y[ index ] = dw_dx * vv;

			    --index;
		  }
	}
}

void fadeR()
{
	const double holdAmount = 1 - _fadeSpeed;

	_avgDensity = 0;
	_avgSpeed = 0;

	double totalDeviations = 0, currentDeviation, tmp_r;
	register double *uop, *vop, *rop, *rp;

	int i = -1;
	uop = uOld;
	vop = vOld;
	rop = rOld;
	rp = r;
	while ( ++i < numCells )
	{
		*(uop++) = *(vop++) = *(rop++) = 0.0;

		_avgSpeed += u[i] * u[i] + v[i] * v[i];

		tmp_r = FMIN(1.0, *(rp++));
		_avgDensity += tmp_r;

		currentDeviation = tmp_r - _avgDensity;
		totalDeviations += currentDeviation * currentDeviation;

		r[i] = tmp_r * holdAmount;
	}
	_avgDensity *= invNumCells;

	_uniformity = 1.0 / (1 + totalDeviations * invNumCells);
}

void fadeRGB()
{
	const double holdAmount = 1 - _fadeSpeed;

	_avgDensity = 0;
	_avgSpeed = 0;

	double totalDeviations = 0.0, currentDeviation, density;
	double tmp_r, tmp_g, tmp_b;

	register double *uop, *vop, *rop, *gop, *bop, *rp, *gp, *bp, *up, *vp;

	uop = uOld;
	vop = vOld;
	rop = rOld;
	gop = gOld;
	bop = bOld;
	gp = g;
	bp = b;
	for ( rp = r, up = u, vp = v; rp < r+numCells; )
	{
		*(uop++) = *(vop++) = *(rop++) = *(gop++) = *(bop++) = 0.0;

		_avgSpeed += *(up) * *(up) + *(vp) * *(vp);

		tmp_r = FMIN(1.0, *(rp));
		tmp_g = FMIN(1.0, *(gp));
		tmp_b = FMIN(1.0, *(bp));

		density = FMAX(tmp_r, FMAX(tmp_g, tmp_b));
		_avgDensity += density;

		currentDeviation = density - _avgDensity;
		totalDeviations += currentDeviation * currentDeviation;

		*(rp++) = tmp_r * holdAmount;
		*(gp++) = tmp_g * holdAmount;
		*(bp++) = tmp_b * holdAmount;

		up++;
		vp++;
	}
	_avgDensity *= invNumCells;
	_avgSpeed *= invNumCells;

	_uniformity = 1.0 / (1 + totalDeviations * invNumCells);
}

void advect(int b, register double *_d, register double *d0, register double *du, register double *dv)
{
	int i, j, i0, j0, i1, j1, index;
	double x, y, s0, t0, s1, t1, dt0x, dt0y;

	dt0x = _dt * NX;
	dt0y = _dt * NY;

	for (j = NY; j > 0; --j)
	{
		index = FLUID_IX(NX, j);
		for (i = NX; i > 0; --i)
		{
			x = i - dt0x * du[index];
			y = j - dt0y * dv[index];

			if (x > NX + 0.5) x = NX + 0.5;
			if (x < 0.5) x = 0.5;

			i0 = (int)x;
			i1 = i0 + 1;

			if (y > NY + 0.5) y = NY + 0.5;
			if (y < 0.5) y = 0.5;

			j0 = (int)y;
			j1 = j0 + 1;

			s1 = x - i0;
			s0 = 1 - s1;
			t1 = y - j0;
			t0 = 1 - t1;

			_d[index] = s0 * (t0 * d0[FLUID_IX(i0, j0)] + t1 * d0[FLUID_IX(i0, j1)]) + s1 * (t0 * d0[FLUID_IX(i1, j0)] + t1 * d0[FLUID_IX(i1, j1)]);
			--index;
		}
	}
	setBoundary(b, _d);
}

void advectRGB(register double *du, register double *dv)
{
	int i, j, i0, j0;
	double x, y, s0, t0, s1, t1, dt0x, dt0y;
	int index;

	dt0x = _dt * NX;
	dt0y = _dt * NY;

	for (j = NY; j > 0; --j)
	{
		index = FLUID_IX(NX, j);
		for (i = NX; i > 0; --i)
		{
			x = i - dt0x * du[index];
			y = j - dt0y * dv[index];

			if (x > NX + 0.5) x = NX + 0.5;
			if (x < 0.5)     x = 0.5;

			i0 = (int)x;

			if (y > NY + 0.5) y = NY + 0.5;
			if (y < 0.5)     y = 0.5;

			j0 = (int)y;

			s1 = x - i0;
			s0 = 1 - s1;
			t1 = y - j0;
			t0 = 1 - t1;

			i0 = FLUID_IX(i0, j0);
            j0 = i0 + NX2;
            r[index] = s0 * ( t0 * rOld[i0] + t1 * rOld[j0] ) + s1 * ( t0 * rOld[i0+1] + t1 * rOld[j0+1] );
            g[index] = s0 * ( t0 * gOld[i0] + t1 * gOld[j0] ) + s1 * ( t0 * gOld[i0+1] + t1 * gOld[j0+1] );
            b[index] = s0 * ( t0 * bOld[i0] + t1 * bOld[j0] ) + s1 * ( t0 * bOld[i0+1] + t1 * bOld[j0+1] );
			--index;
		}
	}
	setBoundaryRGB();
}

void project(register double *x, register double *y, register double *p, register double *div)
{
	int i, j, index;

	double const h = -0.5 / NX;

	for (j = NY; j > 0; --j)
    {
		index = FLUID_IX(NX, j);
		for (i = NX; i > 0; --i)
		{
			div[index] = h * ( x[index+1] - x[index-1] + y[index+NX2] - y[index-NX2] );
			p[index] = 0;
			--index;
		}
    }

	setBoundary(0, div);
	setBoundary(0, p);

	linearSolver(0, p, div, 1.0, 4.0);

	double const fx = 0.5 * NX;
	double const fy = 0.5 * NY;
	for (j = NY; j > 0; --j)
	{
		index = FLUID_IX(NX, j);
		for (i = NX; i > 0; --i)
		{
			x[index] -= fx * (p[index+1] - p[index-1]);
			y[index] -= fy * (p[index+NX2] - p[index-NX2]);
			--index;
		}
	}

	setBoundary(1, x);
	setBoundary(2, y);
}

void linearSolver(int b, register double *x, register double *x0, double a, double c)
{
	int k, i, j, index;

	if( a == 1. && c == 4. )
	{
		for (k = 0; k < _solverIterations; ++k)
		{
			for (j = NY; j > 0 ; --j)
			{
				index = FLUID_IX(NX, j);
				for (i = NX; i > 0 ; --i)
				{
					x[index] = ( x[index-1] + x[index+1] + x[index - NX2] + x[index + NX2] + x0[index] ) * 0.25;
					--index;
				}
			}
			setBoundary( b, x );
		}
	}
	else
	{
		c = 1.0 / c;
		for (k = 0; k < _solverIterations; ++k)
		{
			for (j = NY; j > 0 ; --j)
			{
				index = FLUID_IX(NX, j);
				for (i = NX; i > 0 ; --i)
				{
					x[index] = ( ( x[index-1] + x[index+1] + x[index - NX2] + x[index + NX2] ) * a + x0[index] ) * c;
					--index;
				}
			}
			setBoundary( b, x );
		}
	}
}

void linearSolverRGB(double a, double c)
{
	int k, i, j, index3, index4, index;

	c = 1.0 / c;

	for ( k = 0; k < _solverIterations; ++k )
	{
	    for (j = NY; j > 0; --j)
	    {
	            index = FLUID_IX(NX, j);
				index3 = index - NX2;
				index4 = index + NX2;
				for (i = NX; i > 0; --i)
				{
					r[index] = ( ( r[index-1] + r[index+1]  +  r[index3] + r[index4] ) * a  +  rOld[index] ) * c;
					g[index] = ( ( g[index-1] + g[index+1]  +  g[index3] + g[index4] ) * a  +  gOld[index] ) * c;
					b[index] = ( ( b[index-1] + b[index+1]  +  b[index3] + b[index4] ) * a  +  bOld[index] ) * c;

					--index;
					--index3;
					--index4;
				}
		}
		setBoundaryRGB();
	}
}

void linearSolverUV(double a, double c)
{
	int index, k, i, j;
	c = 1.0 / c;
	for (k = 0; k < _solverIterations; ++k) {
		for (j = NY; j > 0; --j) {
			index = FLUID_IX(NX, j);
			for (i = NX; i > 0; --i) {
				u[index] = ( ( u[index-1] + u[index+1] + u[index - NX2] + u[index + NX2] ) * a  +  uOld[index] ) * c;
				v[index] = ( ( v[index-1] + v[index+1] + v[index - NX2] + v[index + NX2] ) * a  +  vOld[index] ) * c;
				--index;
			}
		}
		setBoundary( 1, u );
        setBoundary( 2, v );
	}
}

void setBoundary(int bound, register double *x)
{
	int dst1, dst2, src1, src2, i;
	const int step = FLUID_IX(0, 1) - FLUID_IX(0, 0);

	dst1 = FLUID_IX(0, 1);
	src1 = FLUID_IX(1, 1);
	dst2 = FLUID_IX(NX+1, 1 );
	src2 = FLUID_IX(NX, 1);

	if( wrap_x == 1 ) {
		src1 ^= src2;
		src2 ^= src1;
		src1 ^= src2;
	}
	if( bound == 1 && wrap_x == 0 ) {
		for (i = NY; i > 0; --i )
		{
			x[dst1] = -x[src1];     dst1 += step;   src1 += step;
			x[dst2] = -x[src2];     dst2 += step;   src2 += step;
		}
	} else {
		for (i = NY; i > 0; --i )
		{
			x[dst1] = x[src1];      dst1 += step;   src1 += step;
			x[dst2] = x[src2];      dst2 += step;   src2 += step;
		}
	}

	dst1 = FLUID_IX(1, 0);
	src1 = FLUID_IX(1, 1);
	dst2 = FLUID_IX(1, NY+1);
	src2 = FLUID_IX(1, NY);

	if( wrap_y == 1 ) {
		src1 ^= src2;
		src2 ^= src1;
		src1 ^= src2;
	}
	if( bound == 2 && wrap_y == 0 ) {
		for (i = NX; i > 0; --i )
		{
		        x[dst1++] = -x[src1++];
		        x[dst2++] = -x[src2++];
		}
	} else {
		for (i = NX; i > 0; --i )
		{
		        x[dst1++] = x[src1++];
		        x[dst2++] = x[src2++];
		}
	}

	x[FLUID_IX(  0,   0)] = 0.5 * (x[FLUID_IX(1, 0  )] + x[FLUID_IX(  0, 1)]);
	x[FLUID_IX(  0, NY+1)] = 0.5 * (x[FLUID_IX(1, NY+1)] + x[FLUID_IX(  0, NY)]);
	x[FLUID_IX(NX+1,   0)] = 0.5 * (x[FLUID_IX(NX, 0  )] + x[FLUID_IX(NX+1, 1)]);
	x[FLUID_IX(NX+1, NY+1)] = 0.5 * (x[FLUID_IX(NX, NY+1)] + x[FLUID_IX(NX+1, NY)]);

}

void setBoundaryRGB()
{
	if( wrap_x == 0 && wrap_y == 0 ) return;

	int dst1, dst2, src1, src2, i;
	const int step = FLUID_IX(0, 1) - FLUID_IX(0, 0);

	if ( wrap_x == 1 ) {
		dst1 = FLUID_IX(0, 1);
		src1 = FLUID_IX(1, 1);
		dst2 = FLUID_IX(NX+1, 1 );
		src2 = FLUID_IX(NX, 1);

		src1 ^= src2;
		src2 ^= src1;
		src1 ^= src2;

		for (i = NY; i > 0; --i )
		{
			r[dst1] = r[src1]; g[dst1] = g[src1]; b[dst1] = b[src1]; dst1 += step;   src1 += step;
			r[dst2] = r[src2]; g[dst2] = g[src2]; b[dst2] = b[src2]; dst2 += step;   src2 += step;
		}
	}

	if ( wrap_y == 1 ) {
		dst1 = FLUID_IX(1, 0);
		src1 = FLUID_IX(1, 1);
		dst2 = FLUID_IX(1, NY+1);
		src2 = FLUID_IX(1, NY);

		src1 ^= src2;
		src2 ^= src1;
		src1 ^= src2;

		for (i = NX; i > 0; --i )
		{
			r[dst1] = r[src1]; g[dst1] = g[src1]; b[dst1] = b[src1];  ++dst1; ++src1;
			r[dst2] = r[src2]; g[dst2] = g[src2]; b[dst2] = b[src2];  ++dst2; ++src2;
		}
	}
}

AS3_Val setupSolver(void* self, AS3_Val args)
{
	AS3_ArrayValue(args, "IntType, IntType, IntType, IntType", &NX, &NY, &screenW, &screenH);

	numCells = (NX + 2) * (NY + 2);
	NX2 = NX + 2;
	NY2 = NY + 2;

	invNX = 1.0 / NX;
	invNY = 1.0 / NY;
	invNumCells = 1.0 / numCells;

	isw = (double)(1.0 / screenW);
	ish = (double)(1.0 / screenH);

	_dt = FLUID_DEFAULT_DT;
	_fadeSpeed = FLUID_DEFAULT_FADESPEED;
	_visc = FLUID_DEFAULT_VISC;
	_solverIterations = FLUID_DEFAULT_SOLVER_ITERATIONS;
	_colorDiffusion = FLUID_DEFAULT_COLOR_DIFFUSION;
	_doVorticityConfinement = FLUID_DEFAULT_VORTICITY_CONFINEMENT;

	reset();

	return AS3_Ptr(fluidsImage);
}

AS3_Val updateSolver(void* self, AS3_Val args)
{
	addSourceUV();

	if( _doVorticityConfinement )
	{
		calcVorticityConfinement(uOld, vOld);
		addSourceUV();
	}

	swapUV();

	diffuseUV(_visc);

	project(u, v, uOld, vOld);

	swapUV();

	advect(1, u, uOld, uOld, vOld);
	advect(2, v, vOld, uOld, vOld);

	project(u, v, uOld, vOld);

	if(_isRGB == 1) {
		addSourceRGB();
		swapRGB();

		if( _colorDiffusion != 0 && _dt != 0 )
        {
			diffuseRGB(_colorDiffusion);
			swapRGB();
        }

		advectRGB(u, v);

		fadeRGB();
	} else {
		addSource(r, rOld);
		swapR();

		if( _colorDiffusion != 0 && _dt != 0 )
        {
			diffuse(0, r, rOld, _colorDiffusion);
			swapRGB();
        }

		advect(0, r, rOld, u, v);
		fadeR();
	}

	register double *pp;
	for(pp = particlesPool; pp < particlesPool+20;)
	{
		double px = *(pp++);
		double py = *(pp++);
		if(px > 0. || py > 0.)
		{
			addParticles(px, py, 20);
		}
		*(pp-1) = *(pp-2) = 0.0;
	}

	if(drawMode == 0) {
		drawFluidImage();
	} else if(drawMode == 1 || drawMode == 2){
		drawFluidImage();
		drawParticleImage();
	} else if(drawMode == 3 || drawMode == 4){
		drawParticleImage();
	}

	return 0;
}

AS3_Val getParticlesCountPointer(void* self, AS3_Val args)
{
	return AS3_Ptr(&particlesNum);
}

AS3_Val clearParticles(void* self, AS3_Val args)
{
	memset(particles, 0.0, PARTICLES_MAX*6*sizeof(double));
	memset(particles2, 0.0, PARTICLES_MAX*6*sizeof(double));
	_particles = particles;
	_particles2 = particles2;
	particlesNum = 0;

	return 0;
}

AS3_Val setDrawMode(void* self, AS3_Val args)
{
	AS3_ArrayValue(args, "IntType", &drawMode);
	return 0;
}

AS3_Val getParticlesPointer(void* self, AS3_Val args)
{
	return AS3_Ptr(particlesImage);
}
AS3_Val getParticlesPoolPointer(void* self, AS3_Val args)
{
	return AS3_Ptr(particlesPool);
}
AS3_Val getParticlesDataPointer(void* self, AS3_Val args)
{
	return AS3_Ptr(&_particles);
}

AS3_Val getROldPointer(void* self, AS3_Val args)
{
	return AS3_Ptr(rOld);
}
AS3_Val getGOldPointer(void* self, AS3_Val args)
{
	return AS3_Ptr(gOld);
}
AS3_Val getBOldPointer(void* self, AS3_Val args)
{
	return AS3_Ptr(bOld);
}
AS3_Val getUOldPointer(void* self, AS3_Val args)
{
	return AS3_Ptr(uOld);
}
AS3_Val getVOldPointer(void* self, AS3_Val args)
{
	return AS3_Ptr(vOld);
}

AS3_Val setWrap(void* self, AS3_Val args)
{
	AS3_ArrayValue(args, "IntType, IntType", &wrap_x, &wrap_y);
	return 0;
}

AS3_Val setcolorDiffusion(void* self, AS3_Val args)
{
	AS3_ArrayValue(args, "DoubleType", &_colorDiffusion);
	return 0;
}

AS3_Val setsolverIterations(void* self, AS3_Val args)
{
	AS3_ArrayValue(args, "IntType", &_solverIterations);
	return 0;
}

AS3_Val setVorticityConfinement(void* self, AS3_Val args)
{
	AS3_ArrayValue(args, "IntType", &_doVorticityConfinement);
	return 0;
}

AS3_Val setFadeSpeed(void* self, AS3_Val args)
{
	AS3_ArrayValue(args, "DoubleType", &_fadeSpeed);
	return 0;
}

AS3_Val setViscosity(void* self, AS3_Val args)
{
	AS3_ArrayValue(args, "DoubleType", &_visc);
	return 0;
}

AS3_Val setDeltaT(void* self, AS3_Val args)
{
	AS3_ArrayValue(args, "DoubleType", &_dt);
	return 0;
}

int main()
{
	AS3_Val setupSolver_m = AS3_Function( NULL, setupSolver );
	AS3_Val updateSolver_m = AS3_Function( NULL, updateSolver );
	AS3_Val getParticlesPointer_m = AS3_Function( NULL, getParticlesPointer );
	AS3_Val getParticlesPoolPointer_m = AS3_Function( NULL, getParticlesPoolPointer );
	AS3_Val getROldPointer_m = AS3_Function( NULL, getROldPointer );
	AS3_Val getGOldPointer_m = AS3_Function( NULL, getGOldPointer );
	AS3_Val getBOldPointer_m = AS3_Function( NULL, getBOldPointer );
	AS3_Val getUOldPointer_m = AS3_Function( NULL, getUOldPointer );
	AS3_Val getVOldPointer_m = AS3_Function( NULL, getVOldPointer );

	AS3_Val setWrap_m = AS3_Function( NULL, setWrap );
	AS3_Val setcolorDiffusion_m = AS3_Function( NULL, setcolorDiffusion );
	AS3_Val setsolverIterations_m = AS3_Function( NULL, setsolverIterations );
	AS3_Val setFadeSpeed_m = AS3_Function( NULL, setFadeSpeed );
	AS3_Val setViscosity_m = AS3_Function( NULL, setViscosity );
	AS3_Val setDeltaT_m = AS3_Function( NULL, setDeltaT );
	AS3_Val setDrawMode_m = AS3_Function( NULL, setDrawMode );
	AS3_Val setVorticityConfinement_m = AS3_Function( NULL, setVorticityConfinement );
	AS3_Val getParticlesCountPointer_m = AS3_Function( NULL, getParticlesCountPointer );
	AS3_Val clearParticles_m = AS3_Function( NULL, clearParticles );
	AS3_Val getParticlesDataPointer_m = AS3_Function( NULL, getParticlesDataPointer );



	AS3_Val result = AS3_Object("setupSolver: AS3ValType, updateSolver: AS3ValType, getParticlesPointer: AS3ValType, getParticlesPoolPointer: AS3ValType, getROldPointer: AS3ValType, getGOldPointer: AS3ValType, getBOldPointer: AS3ValType, getUOldPointer: AS3ValType, getVOldPointer: AS3ValType, setWrap: AS3ValType, setcolorDiffusion: AS3ValType, setsolverIterations: AS3ValType, setFadeSpeed: AS3ValType, setViscosity: AS3ValType, setDeltaT: AS3ValType, setDrawMode: AS3ValType, setVorticityConfinement: AS3ValType, getParticlesCountPointer: AS3ValType, clearParticles: AS3ValType, getParticlesDataPointer: AS3ValType",
								setupSolver_m, updateSolver_m, getParticlesPointer_m, getParticlesPoolPointer_m, getROldPointer_m,
								getGOldPointer_m, getBOldPointer_m, getUOldPointer_m, getVOldPointer_m,
								setWrap_m, setcolorDiffusion_m, setsolverIterations_m, setFadeSpeed_m, setViscosity_m, setDeltaT_m,
								setDrawMode_m, setVorticityConfinement_m, getParticlesCountPointer_m, clearParticles_m,
								getParticlesDataPointer_m);

	AS3_Release( setupSolver_m );
	AS3_Release( updateSolver_m );
	AS3_Release( getParticlesPointer_m );
	AS3_Release( getParticlesPoolPointer_m );
	AS3_Release( getROldPointer_m );
	AS3_Release( getGOldPointer_m );
	AS3_Release( getBOldPointer_m );
	AS3_Release( getUOldPointer_m );
	AS3_Release( getVOldPointer_m );
	AS3_Release( setWrap_m );
	AS3_Release( setcolorDiffusion_m );
	AS3_Release( setsolverIterations_m );
	AS3_Release( setFadeSpeed_m );
	AS3_Release( setViscosity_m );
	AS3_Release( setDeltaT_m );
	AS3_Release( setDrawMode_m );
	AS3_Release( setVorticityConfinement_m );
	AS3_Release( getParticlesCountPointer_m );
	AS3_Release( clearParticles_m );
	AS3_Release( getParticlesDataPointer_m );

	AS3_LibInit( result );

	return 0;
}