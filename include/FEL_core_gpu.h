/*
 * FEL_core_gpu.h
 *
 *  Created on: Jul 27, 2012
 *      Author: chaudhua
 */

#ifndef FEL_CORE_GPU_H_
#define FEL_CORE_GPU_H_

#include <cstdio>

//#include <cutil.h>

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


void kernel_wrapper( double** rangeData, double** domainData, int nRange, int nDomain, int nBin );
void kernel_wrapper_incremental_1( double** domainData, int nDomain, int nBin );
void kernel_wrapper_incremental_2( double** rangeData, int nBin );


#ifdef __CUDACC__
__global__ void encode( double* rangeData_d, double* domainData_d,
						int* matchedDomainIdList_d,
						double* latestErrorList_d,
						int nRange, int nDomain, int nBin,
						int domStartID
						);
__global__ void encode_incremental( double* rangeData_d, double* domainData_d,
									int* matchedDomainIdList_d,
									double* latestErrorList_d,
									int nRange, int nDomain, int nBin,
									int domStartID
								  );
__global__ void decode( float *pDataA, float *pDataB, float *pResult);
#endif


#endif /* FEL_CORE_GPU_H_ */
