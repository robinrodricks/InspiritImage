#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "AS3.h"


float *realData;
float *imagData;
float *realFFTData;
float *imagFFTData;
float *amplFFTData;
float *phaseFFTData;

float *shiftedData;
int *drawData;

int imageW;
int imageH;
int imageW2;
int imageH2;
int numChannels;
int area2;

int spectrumLength;

static const float INV_LOG2 = 1.4426950408889634;

static void cleanUp();
static void calculateAmp(int len, float *gAmp, float *gRe, float *gIm);
static void calculateSpectrum(int len, float *gAmp, float *gRe, float *gIm, int normalize);
static void calculatePhase(int len, float *gPhase, float *gRe, float *gIm);
static void FFT2DRGB(int n, int m, int inverse, float *gRe, float *gIm, float *GRe, float *GIm);
static void FFT2DGray(int n, int m, int inverse, float *gRe, float *gIm, float *GRe, float *GIm);
static void FFT1D(int n, int inverse, float *gRe, float *gIm, float *GRe, float *GIm);
static void plotImageData(int tw, int th, float *data, int *outData, int shift, int add128);
static void shiftData(float *from, float *to);
static void plotImageCarefull(float *data, int *outData, int shift);


int sampleRate = 44100;
int octaves;
int avgPerOctave;
int averagesNumber = 256;
int spectrumAverageType = 0;
float bandWidth;

//static float indexToFreq(int i);
static int freqToIndex(float freq);
static float calcAvg(int offset, float lowFreq, float hiFreq);


static AS3_Val getBufferPointers(void* self, AS3_Val args)
{
	AS3_Val pointers = AS3_Array("AS3ValType", NULL);
	
	AS3_Set(pointers, AS3_Int(0), AS3_Ptr(&realData));
	AS3_Set(pointers, AS3_Int(1), AS3_Ptr(&imagData));
	AS3_Set(pointers, AS3_Int(2), AS3_Ptr(&realFFTData));
	AS3_Set(pointers, AS3_Int(3), AS3_Ptr(&imagFFTData));
	AS3_Set(pointers, AS3_Int(4), AS3_Ptr(&amplFFTData));
	AS3_Set(pointers, AS3_Int(5), AS3_Ptr(&drawData));
	AS3_Set(pointers, AS3_Int(6), AS3_Ptr(&shiftedData));
	AS3_Set(pointers, AS3_Int(7), AS3_Ptr(&phaseFFTData));
	
	return pointers;
}

static AS3_Val allocateBuffers(void* self, AS3_Val args)
{
	int nw, nh, nw2, nh2, numCh;
	AS3_ArrayValue(args, "IntType, IntType, IntType, IntType, IntType", &nw, &nh, &nw2, &nh2, &numCh );
	
	if(nw != imageW || nh != imageH || numCh != numChannels)
	{
		imageW = nw;
		imageH = nh;
		imageW2 = nw2;
		imageH2 = nh2;
		numChannels = numCh;
		
		cleanUp();

		area2 = (imageW2 * imageH2) * numChannels;

		realData = (float*)malloc( area2 * sizeof(float) );
		imagData = (float*)malloc( area2 * sizeof(float) );

		realFFTData = (float*)malloc( area2 * sizeof(float) );
		imagFFTData = (float*)malloc( area2 * sizeof(float) );

		amplFFTData = (float*)malloc( area2 * sizeof(float) );
		phaseFFTData = (float*)malloc( area2 * sizeof(float) );
		
		shiftedData = (float*)malloc( area2 * sizeof(float) );
		
		drawData = (int*)malloc( (imageW2 * imageH2) * sizeof(int) );
	}
	
	memset(realData, 0.0, area2 * sizeof(float));
	memset(imagData, 0.0, area2 * sizeof(float));
	
	return 0;
}

static AS3_Val initSignalBuffers(void* self, AS3_Val args)
{
	int nl, nch;
	AS3_ArrayValue( args, "IntType, IntType", &nl, &nch );
	
	if(nl != imageW || nch != numChannels)
	{
		imageW = imageW2 = nl;
		numChannels = nch;
		
		bandWidth = (2.0 / (float)imageW) * ((float)sampleRate / 2.0);
		spectrumLength = (imageW >> 1) + 1;
		
		cleanUp();

		area2 = imageW * numChannels;
		
		realData = (float*)malloc( area2 * sizeof(float) );
		imagData = (float*)malloc( area2 * sizeof(float) );

		realFFTData = (float*)malloc( area2 * sizeof(float) );
		imagFFTData = (float*)malloc( area2 * sizeof(float) );

		amplFFTData = (float*)malloc( area2 * sizeof(float) );
		phaseFFTData = (float*)malloc( area2 * sizeof(float) );
		
		shiftedData = (float*)malloc( area2 * sizeof(float) );
	}
	
	memset(realData, 0.0, area2 * sizeof(float));
	memset(imagData, 0.0, area2 * sizeof(float));
	
	return 0;
}

static AS3_Val splitChannels(void* self, AS3_Val args)
{
	int offset = imageW;
	register float *re1 = realData, *re2 = realData+offset, *input = shiftedData, *end = realData+offset;
	
	for(; re1 < end;)
	{
		*(re1++) = *(input++);
		*(re2++) = *(input++);
	}
	
	return 0;
}

static AS3_Val mergeChannels(void* self, AS3_Val args)
{
	int offset = imageW;
	register float *re1 = realData, *re2 = realData+offset, *input = shiftedData, *end = realData+offset;
	
	for(; re1 < end;)
	{
		*(input++) = *(re1++);
		*(input++) = *(re2++);
	}
	
	return 0;
}

static AS3_Val initNoAverages(void* self, AS3_Val args)
{	
	spectrumAverageType = 0;
	
	return 0;
}
static AS3_Val initLinearAverages(void* self, AS3_Val args)
{
	AS3_ArrayValue( args, "IntType", &averagesNumber );
	
	spectrumAverageType = 1;
	
	return 0;
}
static AS3_Val initLogarithmicAverages(void* self, AS3_Val args)
{
	int minBandwidth = 11;
	int bandsPerOctave = 1;
	
	AS3_ArrayValue( args, "IntType, IntType, IntType", &minBandwidth, &bandsPerOctave, &sampleRate );
	
	bandWidth = (2.0 / (float)imageW) * ((float)sampleRate / 2.0);
	
	float nyq = (float)sampleRate / 2.0;
	octaves = 1;
	while ((nyq /= 2.0) > minBandwidth)
	{
		octaves++;
	}
	
	avgPerOctave = bandsPerOctave;
	averagesNumber = octaves * bandsPerOctave;
	
	spectrumAverageType = 2;
	
	return AS3_Int(averagesNumber);
}

static AS3_Val analyzeSpectrum(void* self, AS3_Val args)
{
	int normalizeSpectrum = 0;
	
	AS3_ArrayValue( args, "IntType", &normalizeSpectrum );
	
	calculateSpectrum(spectrumLength, amplFFTData, realFFTData, imagFFTData, normalizeSpectrum);
	
	if (spectrumAverageType == 1)
    {
		int avgWidth = (int) spectrumLength / averagesNumber;
		int i, j, offset;
		float avg;
		for (i = 0; i < averagesNumber; i++)
		{
			avg = 0.0;
			for (j = 0; j < avgWidth; j++)
			{
				offset = j + i * avgWidth;
				if (offset < spectrumLength)
				{
					avg += amplFFTData[offset];
				}
				else
				{
					break;
				}
			}
			shiftedData[i] = avg / (float)(j+1);
		}
		
		if(numChannels == 2)
		{
			for (i = 0; i < averagesNumber; i++)
			{
				avg = 0.0;
				for (j = 0; j < avgWidth; j++)
				{
					offset = j + i * avgWidth;
					if (offset < spectrumLength)
					{
						avg += amplFFTData[ spectrumLength + offset ];
					}
					else
					{
						break;
					}
				}
				shiftedData[ averagesNumber + i ] = avg / (float)(j+1);
			}
		}
	}
	else if (spectrumAverageType == 2)
    {
		int i, j, offset;
		float lowFreq, hiFreq, freqStep, f;
		for (i = 0; i < octaves; i++)
		{
		
			if (i == 0)
			{
				lowFreq = 0;
			}
			else
			{
				lowFreq = (float)(sampleRate / 2.0) / (float)pow(2.0, octaves - i);
			}
			hiFreq = (float)(sampleRate / 2.0) / (float)pow(2.0, octaves - i - 1);
			freqStep = (hiFreq - lowFreq) / avgPerOctave;
			f = lowFreq;
			for (j = 0; j < avgPerOctave; j++)
			{
				offset = j + i * avgPerOctave;
				shiftedData[offset] = calcAvg(0, f, f + freqStep);
				f += freqStep;
			}
		}
	
		if(numChannels == 2)
		{
			for (i = 0; i < octaves; i++)
			{
				if (i == 0)
				{
					lowFreq = 0;
				}
				else
				{
					lowFreq = (float)(sampleRate / 2.0) / (float)pow(2.0, octaves - i);
				}
				hiFreq = (float)(sampleRate / 2.0) / (float)pow(2.0, octaves - i - 1);
				freqStep = (hiFreq - lowFreq) / avgPerOctave;
				f = lowFreq;
				for (j = 0; j < avgPerOctave; j++)
				{
					offset = j + i * avgPerOctave;
					shiftedData[averagesNumber + offset] = calcAvg(spectrumLength, f, f + freqStep);
					f += freqStep;
				}
			}
		}
	}
	
	return 0;
}

static AS3_Val analyzeSignal(void* self, AS3_Val args)
{
	int doFFT = 0;
	int doInverse = 0;
	int doSpectrum = 0;
	int normalizeSpectrum = 0;
	
	AS3_ArrayValue(args, "IntType, IntType, IntType, IntType", &doFFT, &doSpectrum, &doInverse, &normalizeSpectrum );
	
	if(numChannels == 2)
	{
		if(doFFT == 1)
		{
			FFT1D(imageW, -1, realData, imagData, realFFTData, imagFFTData);
			FFT1D(imageW, -1, realData+imageW, imagData+imageW, realFFTData+imageW, imagFFTData+imageW);
		}
		
		if(doSpectrum == 1)
		{
			calculateSpectrum(spectrumLength, amplFFTData, realFFTData, imagFFTData, normalizeSpectrum);
		}
		
		if(doInverse == 1)
		{
			FFT1D(imageW, 1, realFFTData, imagFFTData, realData, imagData);
			FFT1D(imageW, 1, realFFTData+imageW, imagFFTData+imageW, realData+imageW, imagData+imageW);
		}
	}
	else
	{
		if(doFFT == 1)
		{
			FFT1D(imageW, -1, realData, imagData, realFFTData, imagFFTData);
		}
		
		if(doSpectrum == 1)
		{
			calculateSpectrum(spectrumLength, amplFFTData, realFFTData, imagFFTData, normalizeSpectrum);
		}
		
		if(doInverse == 1)
		{
			FFT1D(imageW, 1, realFFTData, imagFFTData, realData, imagData);
		}
	}
	
	return 0;
}

static AS3_Val analyzeImage(void* self, AS3_Val args)
{
	int doFFT = 0;
	int doInverse = 0;
	int doAmplitude = 0;
	int doPhase = 0;
	
	AS3_ArrayValue(args, "IntType, IntType, IntType, IntType", &doFFT, &doAmplitude, &doPhase, &doInverse );
	
	if(numChannels == 3)
	{
		if(doFFT == 1)
		{
			FFT2DRGB(imageW2, imageH2, -1, realData, imagData, realFFTData, imagFFTData);
		}
		
		if(doAmplitude == 1)
		{
			calculateAmp(area2, amplFFTData, realFFTData, imagFFTData);
		}
		if(doPhase == 1)
		{
			calculatePhase(area2, phaseFFTData, realFFTData, imagFFTData);
		}
	
		if(doInverse == 1)
		{
			FFT2DRGB(imageW2, imageH2, 1, realFFTData, imagFFTData, realData, imagData);
		}
	}
	else
	{
		if(doFFT == 1)
		{
			FFT2DGray(imageW2, imageH2, -1, realData, imagData, realFFTData, imagFFTData);
		}
		
		if(doAmplitude == 1)
		{
			calculateAmp(area2, amplFFTData, realFFTData, imagFFTData);
		}
		if(doPhase == 1)
		{
			calculatePhase(area2, phaseFFTData, realFFTData, imagFFTData);
		}
	
		if(doInverse == 1)
		{
			FFT2DGray(imageW2, imageH2, 1, realFFTData, imagFFTData, realData, imagData);
		}
	}
	
	return 0;
}

static AS3_Val drawImageData(void* self, AS3_Val args)
{
	int tw = 0;
	int th = 0;
	int shift = 0;
	int add128 = 0;
	int dataType = 0;
	
	AS3_ArrayValue(args, "IntType, IntType, IntType, IntType, IntType", &dataType, &tw, &th, &shift, &add128 );
	
	if(dataType == 0)
	{
		plotImageData(tw, th, realData, drawData, shift, add128);
	} else if(dataType == 1)
	{
		plotImageData(tw, th, imagData, drawData, shift, add128);
	} else if(dataType == 2)
	{
		plotImageData(tw, th, realFFTData, drawData, shift, add128);
	} else if(dataType == 3)
	{
		plotImageData(tw, th, imagFFTData, drawData, shift, add128);
	} else if(dataType == 4)
	{
		plotImageData(tw, th, amplFFTData, drawData, shift, add128);
	}
	
	return 0;
}

static AS3_Val drawImagePreserveData(void* self, AS3_Val args)
{
	int shift = 1;
	int dataType = 0;
	
	AS3_ArrayValue(args, "IntType, IntType", &dataType, &shift );
	
	if(dataType == 0)
	{
		plotImageCarefull(realData, drawData, shift);
	} else if(dataType == 1)
	{
		plotImageCarefull(imagData, drawData, shift);
	} else if(dataType == 2)
	{
		plotImageCarefull(realFFTData, drawData, shift);
	} else if(dataType == 3)
	{
		plotImageCarefull(imagFFTData, drawData, shift);
	} else if(dataType == 4)
	{
		plotImageCarefull(amplFFTData, drawData, shift);
	}
	
	return 0;
}

static AS3_Val shiftImageData(void* self, AS3_Val args)
{
	int dataType = 0;
	int reverse = 0;
	
	AS3_ArrayValue(args, "IntType, IntType", &dataType, &reverse );
	
	if(dataType == 0)
	{
		if(reverse)
		{
			shiftData(shiftedData, realData);
		} else {
			shiftData(realData, shiftedData);
		}
	} else if(dataType == 1)
	{
		if(reverse)
		{
			shiftData(shiftedData, imagData);
		} else {
			shiftData(imagData, shiftedData);
		}
	} else if(dataType == 2)
	{
		if(reverse)
		{
			shiftData(shiftedData, realFFTData);
		} else {
			shiftData(realFFTData, shiftedData);
		}
	} else if(dataType == 3)
	{
		if(reverse)
		{
			shiftData(shiftedData, imagFFTData);
		} else {
			shiftData(imagFFTData, shiftedData);
		}
	} else if(dataType == 4)
	{
		if(reverse)
		{
			shiftData(shiftedData, amplFFTData);
		} else {
			shiftData(amplFFTData, shiftedData);
		}
	}
	
	return 0;
}

static void FFT2DRGB(int n, int m, int inverse, float *gRe, float *gIm, float *GRe, float *GIm)
{
	int x, y, ind1, ind2;
	int i, j, k;
	int i1, l1, l2, l;
	float tx = 0, ty = 0;
	float u1, u2, t1, t2, z, ca, sa, d;
	
	int l2n = log(n) * INV_LOG2 + 0.5;
	int l2m = log(m) * INV_LOG2 + 0.5;
	
	//Erase all history of this array
	memcpy(GRe, gRe, area2 * sizeof(float));
	memcpy(GIm, gIm, area2 * sizeof(float));	
   
	//Bit reversal of each row
	for(y = 0; y < m; y++) //for each row
	{
		j = 0;
		for(i = 0; i < n - 1; i++)
		{
			ind1 = (y * n + i) * 3;
			ind2 = (y * n + j) * 3;
			
			// R
			GRe[ind1] = gRe[ind2];
			GIm[ind1++] = gIm[ind2++];
			// G
			GRe[ind1] = gRe[ind2];
			GIm[ind1++] = gIm[ind2++];
			// B
			GRe[ind1] = gRe[ind2];
			GIm[ind1] = gIm[ind2];
			
			k = n>>1;
			while (k <= j) {j -= k; k>>=1;}
			j += k;
		}
	}
	
	//Bit reversal of each column
	for(x = 0; x < n; x++) //for each column
	{
		j = 0;
		for(i = 0; i < m - 1; i++)
		{
			if(i < j)
			{
				ind1 = (i * n + x) * 3;
				ind2 = (j * n + x) * 3;
				
				// R
				tx = GRe[ind1];
				ty = GIm[ind1];
				GRe[ind1] = GRe[ind2];
				GIm[ind1++] = GIm[ind2];
				GRe[ind2] = tx;
				GIm[ind2++] = ty;
				// G
				tx = GRe[ind1];
				ty = GIm[ind1];
				GRe[ind1] = GRe[ind2];
				GIm[ind1++] = GIm[ind2];
				GRe[ind2] = tx;
				GIm[ind2++] = ty;
				// B
				tx = GRe[ind1];
				ty = GIm[ind1];
				GRe[ind1] = GRe[ind2];
				GIm[ind1] = GIm[ind2];
				GRe[ind2] = tx;
				GIm[ind2] = ty;
			}
			k = m>>1;
			while (k <= j) {j -= k; k>>=1;}
			j += k;
		}
	}
	
	//Calculate the FFT of the columns
	for(x = 0; x < n; x++) //for each column
	{
		//This is the 1D FFT:
		ca = -1.0;
		sa = 0.0;
		l1 = 1;
		l2 = 1;
		for(l=0; l < l2n; l++)
		{
			l1 = l2;
			l2 <<= 1;
			u1 = 1.0;
			u2 = 0.0;
			for(j = 0; j < l1; j++)
			{
				for(i = j; i < n; i += l2)
				{
					i1 = i + l1;
					
					ind1 = (i1 * n + x) * 3;
					ind2 = (i * n + x) * 3;
					
					// R
					t1 = u1 * GRe[ind1] - u2 * GIm[ind1];
					t2 = u1 * GIm[ind1] + u2 * GRe[ind1];
					GRe[ind1] = GRe[ind2] - t1;
					GIm[ind1++] = GIm[ind2] - t2;
					GRe[ind2] += t1;
					GIm[ind2++] += t2;
					// G
					t1 = u1 * GRe[ind1] - u2 * GIm[ind1];
					t2 = u1 * GIm[ind1] + u2 * GRe[ind1];
					GRe[ind1] = GRe[ind2] - t1;
					GIm[ind1++] = GIm[ind2] - t2;
					GRe[ind2] += t1;
					GIm[ind2++] += t2;
					// B
					t1 = u1 * GRe[ind1] - u2 * GIm[ind1];
					t2 = u1 * GIm[ind1] + u2 * GRe[ind1];
					GRe[ind1] = GRe[ind2] - t1;
					GIm[ind1] = GIm[ind2] - t2;
					GRe[ind2] += t1;
					GIm[ind2] += t2;
				}
				z =  u1 * ca - u2 * sa;
				u2 = u1 * sa + u2 * ca;
				u1 = z;
			}
			sa = (float)inverse * sqrt((1.0 - ca) * 0.5);
			ca = sqrt((1.0 + ca) * 0.5);
		}
	}
	
	//Calculate the FFT of the rows
	for(y = 0; y < m; y++) //for each row
	{
		//This is the 1D FFT:
		ca = -1.0;
		sa = 0.0;
		l1= 1;
		l2 = 1;
		for(l = 0; l < l2m; l++)
		{
			l1 = l2;
			l2 <<= 1;
			u1 = 1.0;
			u2 = 0.0;
			for(j = 0; j < l1; j++)
			{
				for(i = j; i < n; i += l2)
				{
					i1 = i + l1;
					
					ind1 = (y * n + i1) * 3;
					ind2 = (y * n + i) * 3;
					
					// R
					t1 = u1 * GRe[ind1] - u2 * GIm[ind1];
					t2 = u1 * GIm[ind1] + u2 * GRe[ind1];
					GRe[ind1] = GRe[ind2] - t1;
					GIm[ind1++] = GIm[ind2] - t2;
					GRe[ind2] += t1;
					GIm[ind2++] += t2;
					// G
					t1 = u1 * GRe[ind1] - u2 * GIm[ind1];
					t2 = u1 * GIm[ind1] + u2 * GRe[ind1];
					GRe[ind1] = GRe[ind2] - t1;
					GIm[ind1++] = GIm[ind2] - t2;
					GRe[ind2] += t1;
					GIm[ind2++] += t2;
					// B
					t1 = u1 * GRe[ind1] - u2 * GIm[ind1];
					t2 = u1 * GIm[ind1] + u2 * GRe[ind1];
					GRe[ind1] = GRe[ind2] - t1;
					GIm[ind1] = GIm[ind2] - t2;
					GRe[ind2] += t1;
					GIm[ind2] += t2;
				}
				z =  u1 * ca - u2 * sa;
				u2 = u1 * sa + u2 * ca;
				u1 = z;
			}
			sa = (float)inverse * sqrt((1.0 - ca) * 0.5);
			ca = sqrt((1.0 + ca) * 0.5);
		}
	}
 
	//if(inverse == 1) d = 1.0 / (float)n; else d = 1.0 / (float)m;
	d = 1.0 / (float)n;
	register float *re2, *im2, *end;
	for( re2 = GRe, im2 = GIm, end = GRe+area2; re2 < end; )
	{
		*(re2++) *= d;
		*(im2++) *= d;
	}
}

static void FFT2DGray(int n, int m, int inverse, float *gRe, float *gIm, float *GRe, float *GIm)
{
	int x, y, ind1, ind2;
	int i, j, k;
	int i1, l1, l2, l;
	float tx = 0, ty = 0;
	float u1, u2, t1, t2, z, ca, sa, d;
	
	int l2n = log(n) * INV_LOG2 + 0.5;
	int l2m = log(m) * INV_LOG2 + 0.5;
	
	//Erase all history of this array
	memcpy(GRe, gRe, area2 * sizeof(float));
	memcpy(GIm, gIm, area2 * sizeof(float));	
   
	//Bit reversal of each row
	for(y = 0; y < m; y++) //for each row
	{
		j = 0;
		for(i = 0; i < n - 1; i++)
		{
			ind1 = (y * n + i);
			ind2 = (y * n + j);
			
			GRe[ind1] = gRe[ind2];
			GIm[ind1] = gIm[ind2];
			
			k = n>>1;
			while (k <= j) {j -= k; k>>=1;}
			j += k;
		}
	}
	
	//Bit reversal of each column
	for(x = 0; x < n; x++) //for each column
	{
		j = 0;
		for(i = 0; i < m - 1; i++)
		{
			if(i < j)
			{
				ind1 = (i * n + x);
				ind2 = (j * n + x);
				
				tx = GRe[ind1];
				ty = GIm[ind1];
				GRe[ind1] = GRe[ind2];
				GIm[ind1] = GIm[ind2];
				GRe[ind2] = tx;
				GIm[ind2] = ty;
			}
			k = m>>1;
			while (k <= j) {j -= k; k>>=1;}
			j += k;
		}
	}
	
	//Calculate the FFT of the columns
	for(x = 0; x < n; x++) //for each column
	{
		//This is the 1D FFT:
		ca = -1.0;
		sa = 0.0;
		l1 = 1;
		l2 = 1;
		for(l=0; l < l2n; l++)
		{
			l1 = l2;
			l2 <<= 1;
			u1 = 1.0;
			u2 = 0.0;
			for(j = 0; j < l1; j++)
			{
				for(i = j; i < n; i += l2)
				{
					i1 = i + l1;
					
					ind1 = (i1 * n + x);
					ind2 = (i * n + x);
					
					// R
					t1 = u1 * GRe[ind1] - u2 * GIm[ind1];
					t2 = u1 * GIm[ind1] + u2 * GRe[ind1];
					GRe[ind1] = GRe[ind2] - t1;
					GIm[ind1] = GIm[ind2] - t2;
					GRe[ind2] += t1;
					GIm[ind2] += t2;
				}
				z =  u1 * ca - u2 * sa;
				u2 = u1 * sa + u2 * ca;
				u1 = z;
			}
			sa = (float)inverse * sqrt((1.0 - ca) * 0.5);
			ca = sqrt((1.0 + ca) * 0.5);
		}
	}
	
	//Calculate the FFT of the rows
	for(y = 0; y < m; y++) //for each row
	{
		//This is the 1D FFT:
		ca = -1.0;
		sa = 0.0;
		l1= 1;
		l2 = 1;
		for(l = 0; l < l2m; l++)
		{
			l1 = l2;
			l2 <<= 1;
			u1 = 1.0;
			u2 = 0.0;
			for(j = 0; j < l1; j++)
			{
				for(i = j; i < n; i += l2)
				{
					i1 = i + l1;
					
					ind1 = (y * n + i1);
					ind2 = (y * n + i);
					
					// R
					t1 = u1 * GRe[ind1] - u2 * GIm[ind1];
					t2 = u1 * GIm[ind1] + u2 * GRe[ind1];
					GRe[ind1] = GRe[ind2] - t1;
					GIm[ind1] = GIm[ind2] - t2;
					GRe[ind2] += t1;
					GIm[ind2] += t2;
				}
				z =  u1 * ca - u2 * sa;
				u2 = u1 * sa + u2 * ca;
				u1 = z;
			}
			sa = (float)inverse * sqrt((1.0 - ca) * 0.5);
			ca = sqrt((1.0 + ca) * 0.5);
		}
	}
 
	//if(inverse == 1) d = 1.0 / (float)n; else d = 1.0 / (float)m;
	d = 1.0 / (float)n;
	register float *re2, *im2, *end;
	for( re2 = GRe, im2 = GIm, end = GRe+area2; re2 < end; )
	{
		*(re2++) *= d;
		*(im2++) *= d;
	}
}

static void FFT1D(int n, int inverse, float *gRe, float *gIm, float *GRe, float *GIm)
{
	int i, j, k;
	int i1, l1, l2, l;
	float tx, ty;
	float u1, u2, t1, t2, z, ca, sa, d;
	
	int l2n = log(n) * INV_LOG2 + 0.5;
	
	//Erase all history of this array
	memcpy(GRe, gRe, n * sizeof(float));
	memcpy(GIm, gIm, n * sizeof(float));
	
	j = 0;
	i1 = n>>1;
	for(i = 0; i < n - 1; i++)
	{
		if(i < j)
		{
			tx = GRe[i];
			ty = GIm[i];
			GRe[i] = GRe[j];
			GIm[i] = GIm[j];
			GRe[j] = tx;
			GIm[j] = ty;
		}
		k = i1;
		while (k <= j) {j -= k; k>>=1;}
		j += k;
	}
	
	//This is the 1D FFT:
	ca = -1.0;
	sa = 0.0;
	l1 = 1;
	l2 = 1;
	for(l=0; l < l2n; l++)
	{
		l1 = l2;
		l2 <<= 1;
		u1 = 1.0;
		u2 = 0.0;
		for(j = 0; j < l1; j++)
		{
			for(i = j; i < n; i += l2)
			{
				i1 = i + l1;
				
				t1 = u1 * GRe[i1] - u2 * GIm[i1];
				t2 = u1 * GIm[i1] + u2 * GRe[i1];
				GRe[i1] = GRe[i] - t1;
				GIm[i1] = GIm[i] - t2;
				GRe[i] += t1;
				GIm[i] += t2;
			}
			z =  u1 * ca - u2 * sa;
			u2 = u1 * sa + u2 * ca;
			u1 = z;
		}
		sa = (float)inverse * sqrt((1.0 - ca) * 0.5);
		ca = sqrt((1.0 + ca) * 0.5);
	}
	
	if(inverse == 1) 
	{
		d = 1.0 / (float)n;
	
		register float *re2, *im2, *end;
		for( re2 = GRe, im2 = GIm, end = GRe+n; re2 < end; )
		{
			*(re2++) *= d;
			*(im2++) *= d;
		}
	}
}

static void calculateAmp(int len, float *gAmp, float *gRe, float *gIm)
{
	register float *re, *im, *am, *end;
	for( re = gRe, im = gIm, am = gAmp, end = gAmp+len; am < end; re++, im++)
	{
		*(am++) = sqrt((*re) * (*re) + (*im) * (*im));
	}
}

static void calculateSpectrum(int len, float *gAmp, float *gRe, float *gIm, int normalize)
{
	register float *re, *im, *am, *end;
	float leftpeak = 0.0, rightPeak = 0.0;
	
	for( re = gRe, im = gIm, am = gAmp, end = gAmp+len; am < end; re++, im++, am++)
	{
		*(am) = sqrt((*re) * (*re) + (*im) * (*im));
		if(*am > leftpeak) leftpeak = *am;
	}
	
	if(normalize == 1)
	{
		leftpeak = 1.0 / leftpeak;
		for( am = gAmp, end = gAmp+len; am < end;)
		{
			*(am++) *= leftpeak;
		}
	}
	
	if(numChannels == 2)
	{
		
		for( re = gRe+imageW, im = gIm+imageW, am = gAmp+len, end = gAmp+len+len; am < end; re++, im++, am++)
		{
			*(am) = sqrt((*re) * (*re) + (*im) * (*im));
			if(*am > rightPeak) rightPeak = *am;
		}
		
		if(normalize == 1)
		{
			rightPeak = 1.0 / rightPeak;
			for( am = gAmp+len, end = gAmp+len+len; am < end;)
			{
				*(am++) *= rightPeak;
			}
		}
	}
}

static float calcAvg(int offset, float lowFreq, float hiFreq)
{
	int lowBound = freqToIndex(lowFreq);
	int hiBound = freqToIndex(hiFreq);
	float avg = 0.0;
	int i;
	for (i = lowBound; i <= hiBound; i++)
	{
		avg += amplFFTData[i + offset];
	}
	avg /= (float)(hiBound - lowBound + 1);
	return avg;
}

static int freqToIndex(float freq)
{
	if (freq < bandWidth / 2.0) return 0;
	if (freq > sampleRate / 2.0 - bandWidth / 2.0) return spectrumLength - 1;
	
	float fraction = freq / (float) sampleRate;
	int i = (int)(imageW * fraction + 0.5);
	return i;
}
/*
static float indexToFreq(int i)
{
	float bw = bandWidth;
	
	if ( i == 0 ) return bw * 0.25;
	
	if ( i == spectrumLength - 1 ) 
	{
	  float lastBinBeginFreq = (sampleRate / 2.0) - (bw / 2.0);
	  float binHalfWidth = bw * 0.25;
	  return lastBinBeginFreq + binHalfWidth;
	}
	return i * bw;
}
*/
static void calculatePhase(int len, float *gPhase, float *gRe, float *gIm)
{
	register float *re, *im, *ph, *end;
	for( re = gRe, im = gIm, ph = gPhase, end = gPhase+len; ph < end; re++, im++)
	{
		*(ph++) = atan2((*re), (*im));
	}
}

static void shiftData(float *from, float *to)
{
	int i, j, ind;
	
	int hw = imageW2 >> 1;
	int hh = imageH2 >> 1;
	
	register float *out = to;
	
	if(numChannels == 3)
	{
		for(i = 0; i < imageH2; i++)
		{
			for(j = 0; j < imageW2; j++)
			{
				ind = (((i + hh) % imageH2) * imageW2 + ((j + hw) % imageW2)) * 3;

				*(out++) = from[ind++];
				*(out++) = from[ind++];
				*(out++) = from[ind];
			}
		}
	}
	else
	{
		for(i = 0; i < imageH2; i++)
		{
			for(j = 0; j < imageW2; j++)
			{
				ind = (((i + hh) % imageH2) * imageW2 + ((j + hw) % imageW2));

				*(out++) = from[ind];
			}
		}
	}
}
static void plotImageCarefull(float *data, int *outData, int shift)
{
	int i, j, ind = 0;
	int x, y;
	int r, g, b;
	int hw = imageW2 >> 1;
	int hh = imageH2 >> 1;
	//int stride = (imageW2 - tw) * numChannels;
	
	//float maxd = (float)sqrt( 2.0 * hw * hh );
	float maxd = (float)sqrt( 2.0 * (float)pow( (hw > hh ? hw : hh), 2 ) );
	float d;
	
	register int *out = outData;
	
	if(numChannels == 3)
	{
		for(i = 0; i < imageH2; i++)
		{
			for(j = 0; j < imageW2; j++)
			{
				//if(shift)
				//{
					ind = (((i + hh) % imageH2) * imageW2 + ((j + hw) % imageW2)) * 3;
				//}
				
				x = j-hw;
				y = i-hh;
				
				d = maxd / ( 1.0 +  (float)sqrt( x*x+y*y ) );
				
				r = 127 + ( data[ind] < 0 ? -1 : 1 ) * (int)( abs(data[ind]) / ( 1.7 * d) + 0.5 );
				ind++;
				g = 127 + ( data[ind] < 0 ? -1 : 1 ) * (int)( abs(data[ind]) / ( 1.7 * d) + 0.5 );
				ind++;
				b = 127 + ( data[ind] < 0 ? -1 : 1 ) * (int)( abs(data[ind]) / ( 1.7 * d) + 0.5 );
				ind++;

				r = r < 0 ? 0 : (r > 255 ? 255 : r);
				g = g < 0 ? 0 : (g > 255 ? 255 : g);
				b = b < 0 ? 0 : (b > 255 ? 255 : b);

				*(out++) = ((r<<16)+(g<<8)+b);
			}
			//ind += stride;
		}
	}
	else
	{
		for(i = 0; i < imageH2; i++)
		{
			for(j = 0; j < imageW2; j++)
			{
				//if(shift)
				//{
					ind = (((i + hh) % imageH2) * imageW2 + ((j + hw) % imageW2));
				//}

				x = j-hw;
				y = i-hh;
				
				d = maxd / ( 1.0 +  (float)sqrt( x*x+y*y ) );
				
				r = 127 + ( data[ind] < 0 ? -1.0 : 1.0 ) * (int)( abs(data[ind]) / ( 1.7 * d) + 0.5 );
				ind++;
				
				r = r < 0 ? 0 : (r > 255 ? 255 : r);

				*(out++) = ((r<<16)+(r<<8)+r);
			}
			//ind += stride;
		}
	}
}
static void plotImageData(int tw, int th, float *data, int *outData, int shift, int add128)
{
	int i, j, ind = 0;
	
	int r, g, b;
	int hw = imageW2 >> 1;
	int hh = imageH2 >> 1;
	int stride = (imageW2 - tw) * numChannels;
	
	register int *out = outData;
	
	if(numChannels == 3)
	{
		for(i = 0; i < th; i++)
		{
			for(j = 0; j < tw; j++)
			{
				if(shift)
				{
					ind = (((i + hh) % imageH2) * imageW2 + ((j + hw) % imageW2)) * 3;
				}

				r = (int)data[ind++];
				g = (int)data[ind++];
				b = (int)data[ind++];

				if(add128)
				{
					r += 128;
					g += 128;
					b += 128;
				}

				r = r < 0 ? 0 : (r > 255 ? 255 : r);
				g = g < 0 ? 0 : (g > 255 ? 255 : g);
				b = b < 0 ? 0 : (b > 255 ? 255 : b);

				*(out++) = ((r<<16)+(g<<8)+b);
			}
			ind += stride;
		}
	}
	else
	{
		for(i = 0; i < th; i++)
		{
			for(j = 0; j < tw; j++)
			{
				if(shift)
				{
					ind = (((i + hh) % imageH2) * imageW2 + ((j + hw) % imageW2));
				}

				r = (int)data[ind++];

				if(add128)
				{
					r += 128;
				}

				r = r < 0 ? 0 : (r > 255 ? 255 : r);

				*(out++) = ((r<<16)+(r<<8)+r);
			}
			ind += stride;
		}
	}
}

static void cleanUp()
{
	if(realData) free(realData);
	if(imagData) free(imagData);
	if(realFFTData) free(realFFTData);
	if(imagFFTData) free(imagFFTData);
	if(amplFFTData) free(amplFFTData);
	if(drawData) free(drawData);
	if(shiftedData) free(shiftedData);
	if(phaseFFTData) free(phaseFFTData);
}
static AS3_Val freeBuffers(void* self, AS3_Val args)
{
	cleanUp();
	
	imageW = 0;
	imageH = 0;
	imageW2 = 0;
	imageH2 = 0;
	numChannels = 0;
	
	return 0;
}

int main()
{
	AS3_Val getBufferPointers_m = AS3_Function( NULL, getBufferPointers );
	AS3_Val allocateBuffers_m = AS3_Function( NULL, allocateBuffers );
	AS3_Val freeBuffers_m = AS3_Function( NULL, freeBuffers );
	AS3_Val analyzeImage_m = AS3_Function( NULL, analyzeImage );
	AS3_Val drawImageData_m = AS3_Function( NULL, drawImageData );
	AS3_Val shiftImageData_m = AS3_Function( NULL, shiftImageData );
	AS3_Val drawImagePreserveData_m = AS3_Function( NULL, drawImagePreserveData );
	
	AS3_Val analyzeSignal_m = AS3_Function( NULL, analyzeSignal );
	AS3_Val splitChannels_m = AS3_Function( NULL, splitChannels );
	AS3_Val mergeChannels_m = AS3_Function( NULL, mergeChannels );
	AS3_Val initSignalBuffers_m = AS3_Function( NULL, initSignalBuffers );	
	
	AS3_Val initLinearAverages_m = AS3_Function( NULL, initLinearAverages );
	AS3_Val initLogarithmicAverages_m = AS3_Function( NULL, initLogarithmicAverages );
	AS3_Val initNoAverages_m = AS3_Function( NULL, initNoAverages );
	AS3_Val analyzeSpectrum_m = AS3_Function( NULL, analyzeSpectrum );

	AS3_Val result = AS3_Object("getBufferPointers: AS3ValType, allocateBuffers: AS3ValType, freeBuffers: AS3ValType, analyzeImage: AS3ValType, drawImageData: AS3ValType, shiftImageData: AS3ValType, drawImagePreserveData: AS3ValType, analyzeSignal: AS3ValType, splitChannels: AS3ValType, mergeChannels: AS3ValType, initSignalBuffers: AS3ValType, initLinearAverages: AS3ValType, initLogarithmicAverages: AS3ValType, initNoAverages: AS3ValType, analyzeSpectrum: AS3ValType",
	getBufferPointers_m, allocateBuffers_m, freeBuffers_m, analyzeImage_m, drawImageData_m, shiftImageData_m, drawImagePreserveData_m, analyzeSignal_m, splitChannels_m, mergeChannels_m, initSignalBuffers_m, initLinearAverages_m, initLogarithmicAverages_m, initNoAverages_m, analyzeSpectrum_m);

	AS3_Release( getBufferPointers_m );
	AS3_Release( allocateBuffers_m );
	AS3_Release( freeBuffers_m );
	AS3_Release( analyzeImage_m );
	AS3_Release( drawImageData_m );
	AS3_Release( shiftImageData_m );
	AS3_Release( drawImagePreserveData_m );
	
	AS3_Release( analyzeSignal_m );
	AS3_Release( splitChannels_m );
	AS3_Release( mergeChannels_m );
	AS3_Release( initSignalBuffers_m );
	
	AS3_Release( initLinearAverages_m );
	AS3_Release( initLogarithmicAverages_m );
	AS3_Release( initNoAverages_m );
	AS3_Release( analyzeSpectrum_m );

	AS3_LibInit( result );

	return 0;
}