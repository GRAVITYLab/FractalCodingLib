/*
 * FEL_decoder_gpu.h
 *
 *  Created on: Jun 24, 2013
 *      Author: abon
 */

#ifndef FEL_DECODER_GPU_H_
#define FEL_DECODER_GPU_H_

// CUDA Runtime, Interop, and includes
#include <cuda_runtime.h>
#include <cuda_gl_interop.h>
#include <vector_types.h>
#include <vector_functions.h>
#include <driver_functions.h>

#include <cuda.h>
#include <builtin_types.h>
#include <helper_cuda.h>
//#include <helper_cuda_gl.h>
#include <helper_math.h>
#include <helper_functions.h>
#include <helper_timer.h>

typedef unsigned int  uint;
typedef unsigned char uchar;

typedef struct {
	int4 low;
	int4 high;
} Span;



//typedef unsigned char VolumeType;
typedef unsigned int VolumeType;	//NOTE: here we will store index number of histogram data into the volume, index number is unsigned int
//typedef unsigned short VolumeType;
typedef float HistogramType;
typedef int4 CodebookType;
typedef float TemplatesType;	// NOTE: CUDA kernel doesn't support double (instead, it only support double1, double2), or it will report error when creating channel descriptor
typedef int CodebookTemplateIdType;
typedef float2 ErrorsbookType;

typedef int4 SpanType;
typedef int4 FlexibleCodebookType;
typedef float2 FlexibleErrorsbookType;
typedef int SimpleCountType;
typedef float2 SimpleHistogramType;
typedef float4 FlexTexType;

// flexible block size



//const int nBlocks = (int) 50 * 50 * 10;
//const int nBins = 32;

const int nDimension = 64;
const int flexNBin = 64;
const int nTemplate = 469;	// TODO: WARNING! this number needs to be changed!
const cudaExtent flexibleVolumeSize = make_cudaExtent(64, 64, 32); 	// NOTE: here we just hard code the number of fractal histograms or simple histograms, each is 131072
const cudaExtent flexibleHistogramSize = make_cudaExtent(64, 2048, 64);	// for simpleHistogram, 64*64*32 blocks, 64 bins, width is the bin number
const cudaExtent flexibleTemplatesSize = make_cudaExtent(64, 469, 1); // NOTE: here we hard code the size
//const cudaExtent flexTexSize = make_cudaExtent(3, 3, 3);

//__constant__ float3x4 c_invViewMatrix;  // inverse view matrix

// NOTE: here we assume the MAX number of blocks is 64*64*64=262144.
const int nMaxBlock = 262144;

//__device__ float gsum[4];
//__device__ float td[3][3][3];

__device__ float4 flexBlockData[nMaxBlock];	// x is mean, y is variance, z is entropy

__device__ Span flexBlock[flexNBin * flexNBin * flexNBin]; // span of each block, here we assume the max number of blocks, which means each block has size 1
__device__ int nFlexBlock;		// total number of blocks in this volume
__device__ int nFlexBlockX;		// number of blocks in x dimension
__device__ int nFlexBlockY;
__device__ int nFlexBlockZ;

__device__ int4 corner[8];
__device__ Span subSpan[8][6*6*6];	// 8 corners, each corner has at most 216 subSpans.
__device__ int nSubSpan[8];		// number of sub spans in that corner
__device__ float cornerHistogram[8][6*6*6][flexNBin];
__device__ float cornerSumHistogram[8][flexNBin];	// the sum histogram of all spans of each corner



#ifdef __CUDACC__
// flexible block size
cudaArray *d_codebookSpanLowArray = 0;
cudaArray *d_codebookSpanHighArray = 0;
cudaArray *d_flexibleCodebookArray = 0;
cudaArray *d_flexibleErrorsbookArray = 0;
cudaArray *d_simpleSpanLowArray = 0;
cudaArray *d_simpleSpanHighArray = 0;
cudaArray *d_simpleCountArray = 0;
cudaArray *d_simpleHistogramArray = 0;
cudaArray *d_flexibleTemplatesArray = 0;
cudaArray *d_flexTexArray = 0;

cudaTextureObject_t flexBlockTexObj;
cudaSurfaceObject_t flexBlockSurfObj = 0;

texture<SpanType, cudaTextureType3D, cudaReadModeElementType> codebookSpanLowTex;
texture<SpanType, cudaTextureType3D, cudaReadModeElementType> codebookSpanHighTex;
texture<FlexibleCodebookType, cudaTextureType3D, cudaReadModeElementType> flexibleCodebookTex;
texture<FlexibleErrorsbookType, cudaTextureType2DLayered, cudaReadModeElementType> flexibleErrorsbookTex;
texture<SpanType, cudaTextureType3D, cudaReadModeElementType> simpleSpanLowTex;
texture<SpanType, cudaTextureType3D, cudaReadModeElementType> simpleSpanHighTex;
texture<SimpleCountType, cudaTextureType3D, cudaReadModeElementType> simpleCountTex;
texture<SimpleHistogramType, cudaTextureType2DLayered, cudaReadModeElementType> simpleHistogramTex;
texture<TemplatesType, cudaTextureType2DLayered, cudaReadModeElementType> flexibleTemplatesTex;

texture<FlexTexType, cudaTextureType3D, cudaReadModeElementType> flexBlockTex;

__device__ void flexibleFractalDecoding( float* original, float* decoded, int nBin, int flip, int shift);

__global__ void d_divideBlock(int divPar);
__global__ void d_queryBlockNew(int blockNumber);
__global__ void d_querySpanNew();
__global__ void d_sumSanHistogram(int n);
__global__ void d_computeBlock(int n);
__global__ void d_clearCornerHistogram();
#endif

void initCuda( int4 *h_codebookSpanLow, int4 *h_codebookSpanHigh,
				 int4 *h_flexibleCodebook,
				 float2 *h_flexibleErrorsbook,
				 int4 *h_simpleLow, int4 *h_simpleHigh, int *h_simpleCount,
				 float2 *h_simpleHistogram, float *h_flexibleTemplates );

void dataProcessing();
void freeCudaBuffers();

#endif
/* FEL_DECODER_GPU_H_ */
