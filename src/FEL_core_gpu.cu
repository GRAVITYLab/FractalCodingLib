#include "FEL_core_gpu.h"

void
kernel_wrapper( double** rangeData, double** domainData, int nRange, int nDomain, int nBin )
{
	int* matchedDomainIdList = new int[nRange];
	int* matchedDomainIdList_d = 0;
	double* rangeData_d = 0;
	double* domainData_d = 0;
	double* latestErrorList_d = 0;
	memset( matchedDomainIdList, 0, sizeof(int)*nRange );

	// Initialize threads and blocks
	dim3 threads( 16, 1 );
	dim3 blocks( 16, 12 );

	// Allocate device memory
	//CUDA_SAFE_CALL( 
			cudaMalloc( (void **)&rangeData_d, nRange*nBin*sizeof(double) );// );
	//CUDA_SAFE_CALL( 
			cudaMalloc( (void **)&domainData_d, nDomain*nBin*sizeof(double) );// );
	//CUDA_SAFE_CALL( 
			cudaMalloc( (void **)&matchedDomainIdList_d, nRange*sizeof(int) );// );
	//CUDA_SAFE_CALL( 
			cudaMalloc( (void **)&latestErrorList_d, nRange*sizeof(double) );// );

	// Copy data to device memory
	//CUDA_SAFE_CALL( 
			cudaMemcpy( rangeData_d, (*rangeData), nRange*nBin*sizeof(double), cudaMemcpyHostToDevice );// );
	//CUDA_SAFE_CALL( 
			cudaMemcpy( domainData_d, (*domainData), nDomain*nBin*sizeof(double), cudaMemcpyHostToDevice );// );

	// Launch encoding kernel
	int nDomPerRound = 80;
	int nRound = nDomain / nDomPerRound;
	fprintf( stderr, "Number of rounds needed: %d\n", nRound );

	double gpuTime;
	unsigned int hTimer;
	//CUT_SAFE_CALL( 
			//cutCreateTimer(&hTimer) ;//);
	//CUT_SAFE_CALL( 
			//cutResetTimer(hTimer) ;//);
	//CUT_SAFE_CALL( 
			//cutStartTimer(hTimer) ;//);

	for( int iR = 0; iR<nRound; iR ++ )
	{
		fprintf( stderr, "Kernel call round: %d\n", iR );
		encode<<<blocks, 1>>>( rangeData_d, domainData_d,
							   matchedDomainIdList_d,
							   latestErrorList_d,
							   nRange, nDomain, nBin,
							   iR*nDomPerRound );
	}

	//CUT_CHECK_ERROR("encode() execution failed\n");
	//CUDA_SAFE_CALL( 
		cudaThreadSynchronize() ;//);

	//CUT_SAFE_CALL( 
			//cutStopTimer(hTimer) ;//);
	//gpuTime = cutGetTimerValue(hTimer);

	// Copy data back to device memory
	//CUDA_SAFE_CALL( 
		cudaMemcpy( matchedDomainIdList, matchedDomainIdList_d, nRange*sizeof(int), cudaMemcpyDeviceToHost ) ;//);

	// Deallocate device memory
	//CUDA_SAFE_CALL( 
			cudaFree( latestErrorList_d ) ;//);
	//CUDA_SAFE_CALL( 
			cudaFree( matchedDomainIdList_d ) ;//);
	//CUDA_SAFE_CALL( 
			cudaFree( rangeData_d ) ;//);
	//CUDA_SAFE_CALL( 
			cudaFree( domainData_d ) ;//);

	fprintf( stderr, "Time: %g milliseconds\n", gpuTime );

	//for( int i=0; i<nRange; i++ )
	//	fprintf( stderr, "%d\n", matchedDomainIdList[i] );

	delete [] matchedDomainIdList;

}// end function

/*
void
kernel_wrapper_incremental_1( double** domainData, int nDomain, int nBin )
{
	double* domainData_d = 0;
	
	// Allocate device memory
	CUDA_SAFE_CALL( cudaMalloc( (void **)&domainData_d, nDomain*nBin*sizeof(double) ) );
	
}// end function

void
kernel_wrapper_incremental_2( double** rangeData, int nBin )
{
	int* matchedDomainIdList = new int[nRange];
	int* matchedDomainIdList_d = 0;
	double* rangeData_d = 0;

	double* latestErrorList_d = 0;
	memset( matchedDomainIdList, 0, sizeof(int)*nRange );

	// Initialize threads and blocks
	dim3 threads( 16, 1 );
	dim3 blocks( 16, 12 );

	// Allocate device memory
	CUDA_SAFE_CALL( cudaMalloc( (void **)&rangeData_d, nRange*nBin*sizeof(double) ) );
	CUDA_SAFE_CALL( cudaMalloc( (void **)&matchedDomainIdList_d, nRange*sizeof(int) ) );
	CUDA_SAFE_CALL( cudaMalloc( (void **)&latestErrorList_d, nRange*sizeof(double) ) );

	// Copy data to device memory
	CUDA_SAFE_CALL( cudaMemcpy( rangeData_d, (*rangeData), nRange*nBin*sizeof(double), cudaMemcpyHostToDevice ) );
	CUDA_SAFE_CALL( cudaMemcpy( domainData_d, (*domainData), nDomain*nBin*sizeof(double), cudaMemcpyHostToDevice ) );

	// Launch encoding kernel
	int nDomPerRound = 80;
	int nRound = nDomain / nDomPerRound;
	fprintf( stderr, "Number of rounds needed: %d\n", nRound );

	double gpuTime;
	unsigned int hTimer;
	CUT_SAFE_CALL( cutCreateTimer(&hTimer) );
	CUT_SAFE_CALL( cutResetTimer(hTimer) );
	CUT_SAFE_CALL( cutStartTimer(hTimer) );

	for( int iR = 0; iR<nRound; iR ++ )
	{
		fprintf( stderr, "Kernel call round: %d\n", iR );
		encode<<<blocks, 1>>>( rangeData_d, domainData_d,
							   matchedDomainIdList_d,
							   latestErrorList_d,
							   nRange, nDomain, nBin,
							   iR*nDomPerRound );
	}

	CUT_CHECK_ERROR("encode() execution failed\n");
	CUDA_SAFE_CALL( cudaThreadSynchronize() );

	CUT_SAFE_CALL( cutStopTimer(hTimer) );
	gpuTime = cutGetTimerValue(hTimer);

	// Copy data back to device memory
	CUDA_SAFE_CALL( cudaMemcpy( matchedDomainIdList, matchedDomainIdList_d, nRange*sizeof(int), cudaMemcpyDeviceToHost ) );

	// Deallocate device memory
	CUDA_SAFE_CALL( cudaFree( latestErrorList_d ) );
	CUDA_SAFE_CALL( cudaFree( matchedDomainIdList_d ) );
	CUDA_SAFE_CALL( cudaFree( rangeData_d ) );
	CUDA_SAFE_CALL( cudaFree( domainData_d ) );

	fprintf( stderr, "Time: %g milliseconds\n", gpuTime );

	//for( int i=0; i<nRange; i++ )
	//	fprintf( stderr, "%d\n", matchedDomainIdList[i] );

	delete [] matchedDomainIdList;

}// end function
*/

__global__ void
encode( double* rangeData_d, double* domainData_d,
		int* matchedDomainIdList_d,
		double* latestErrorList_d,
		int nRange, int nDomain, int nBin,
		int domStartID )
{
	int tid = (blockIdx.y * gridDim.x * gridDim.y) + blockIdx.x * gridDim.x + threadIdx.x;

	int bestMatchDomainId;
	double buffer, matchErrorPrcnt, minErrorPrcnt;
	double curRange[130];
	double domainFreqList[130];
	double domainFreqListCopy[130];
	double domainFreqListRev[130];
	double* buffer2[130];

	if( tid < nRange )
	{
		memcpy( curRange, rangeData_d + tid*nBin, sizeof(double)*nBin );

		if( domStartID == 0 )
			minErrorPrcnt = 1000;
		else
			minErrorPrcnt = latestErrorList_d[tid];
		bestMatchDomainId = 0;

		int domMaxId = min ( domStartID+10, nDomain );
		for( int iD = domStartID; iD<domMaxId; iD++ )
		{
			memcpy( domainFreqList, domainData_d + iD*nBin, sizeof(double)*nBin );

			// Circular shift
			memcpy( domainFreqListCopy, domainFreqList, sizeof(double)*nBin );
			for( int iT = 0; iT<nBin; iT++ )
			{
				// Place rightmost to buffer
				buffer = domainFreqListCopy[nBin-1];

				// Move everybody else by 1 bit
				memcpy( buffer2, domainFreqListCopy, sizeof(double)*(nBin-1) );
				memcpy( domainFreqListCopy+1, buffer2, sizeof(double)*(nBin-1) );

				// Fill back the leftmost
				domainFreqListCopy[0] = buffer;

				// Compare with range distribution
				matchErrorPrcnt = 0;
				for( int iB = 0; iB<nBin; iB++ )
				{
					double del = ( domainFreqListCopy[iB] - curRange[iB] );
					matchErrorPrcnt += del * del;
				}

				if( matchErrorPrcnt < minErrorPrcnt )
				{
					minErrorPrcnt = matchErrorPrcnt;
					bestMatchDomainId = iD;
					//optRot = iT;
					//optRef = false;
				}
			}

			// Reflect distribution
			for( int i=0; i<nBin; i++ )
				domainFreqListRev[i] = domainFreqList[nBin-i-1];

			// Circular shift
			memcpy( domainFreqListCopy, domainFreqListRev, sizeof(double)*nBin );
			for( int iT = 0; iT<nBin; iT++ )
			{
				// Place rightmost to buffer
				buffer = domainFreqListCopy[nBin-1];

				// Move everybody else by 1 bit
				memcpy( buffer2, domainFreqListCopy, sizeof(double)*(nBin-1) );
				memcpy( domainFreqListCopy+1, buffer2, sizeof(double)*(nBin-1) );

				// Fill back the leftmost
				domainFreqListCopy[0] = buffer;


				// Compare with range distribution
				matchErrorPrcnt = 0;
				for( int iB = 0; iB<nBin; iB++ )
				{
					double del = ( domainFreqListCopy[iB] - curRange[iB] );
					matchErrorPrcnt += del * del;
				}

				if( matchErrorPrcnt < minErrorPrcnt )
				{
					minErrorPrcnt = matchErrorPrcnt;
					bestMatchDomainId = iD;
					//optRot = iT;
					//optRef = true;
				}
			}
			//(*rotateAmount) = optRot;
			//(*isReflected) = optRef;

		}// end for : scan domains

		matchedDomainIdList_d[tid] = bestMatchDomainId;
		latestErrorList_d[tid] = minErrorPrcnt;

	}// end if (tid < nRange)

}// end function

/*
__global__ void
encode_incremental( double* rangeData_d, double* domainData_d,
					int* matchedDomainIdList_d,
					double* latestErrorList_d,
					int nRange, int nDomain, int nBin,
					int domStartID )
{
	// Declare shared memory which stores 
	// comparison for one range
	double* errorList = new double[nDomain];	
	
	// Each thread can perform one comparison
	int tid = blockIdx.y * gridDim.x * gridDim.y) + blockIdx.x * gridDim.x + threadIdx.x;

	int bestMatchDomainId;
	double buffer, matchErrorPrcnt, minErrorPrcnt;
	double curRange[130];
	double domainFreqList[130];
	double domainFreqListCopy[130];
	double domainFreqListRev[130];
	double* buffer2[130];

	if( tid < nRange )
	{
		memcpy( curRange, rangeData_d + tid*nBin, sizeof(double)*nBin );

		if( domStartID == 0 )
			minErrorPrcnt = 1000;
		else
			minErrorPrcnt = latestErrorList_d[tid];
		bestMatchDomainId = 0;

		int domMaxId = min ( domStartID+10, nDomain );
		for( int iD = domStartID; iD<domMaxId; iD++ )
		{
			memcpy( domainFreqList, domainData_d + iD*nBin, sizeof(double)*nBin );

			// Circular shift
			memcpy( domainFreqListCopy, domainFreqList, sizeof(double)*nBin );
			for( int iT = 0; iT<nBin; iT++ )
			{
				// Place rightmost to buffer
				buffer = domainFreqListCopy[nBin-1];

				// Move everybody else by 1 bit
				memcpy( buffer2, domainFreqListCopy, sizeof(double)*(nBin-1) );
				memcpy( domainFreqListCopy+1, buffer2, sizeof(double)*(nBin-1) );

				// Fill back the leftmost
				domainFreqListCopy[0] = buffer;

				// Compare with range distribution
				matchErrorPrcnt = 0;
				for( int iB = 0; iB<nBin; iB++ )
				{
					double del = ( domainFreqListCopy[iB] - curRange[iB] );
					matchErrorPrcnt += del * del;
				}

				if( matchErrorPrcnt < minErrorPrcnt )
				{
					minErrorPrcnt = matchErrorPrcnt;
					bestMatchDomainId = iD;
					//optRot = iT;
					//optRef = false;
				}
			}

			// Reflect distribution
			for( int i=0; i<nBin; i++ )
				domainFreqListRev[i] = domainFreqList[nBin-i-1];

			// Circular shift
			memcpy( domainFreqListCopy, domainFreqListRev, sizeof(double)*nBin );
			for( int iT = 0; iT<nBin; iT++ )
			{
				// Place rightmost to buffer
				buffer = domainFreqListCopy[nBin-1];

				// Move everybody else by 1 bit
				memcpy( buffer2, domainFreqListCopy, sizeof(double)*(nBin-1) );
				memcpy( domainFreqListCopy+1, buffer2, sizeof(double)*(nBin-1) );

				// Fill back the leftmost
				domainFreqListCopy[0] = buffer;


				// Compare with range distribution
				matchErrorPrcnt = 0;
				for( int iB = 0; iB<nBin; iB++ )
				{
					double del = ( domainFreqListCopy[iB] - curRange[iB] );
					matchErrorPrcnt += del * del;
				}

				if( matchErrorPrcnt < minErrorPrcnt )
				{
					minErrorPrcnt = matchErrorPrcnt;
					bestMatchDomainId = iD;
					//optRot = iT;
					//optRef = true;
				}
			}
			//(*rotateAmount) = optRot;
			//(*isReflected) = optRef;

		}// end for : scan domains

		matchedDomainIdList_d[tid] = bestMatchDomainId;
		latestErrorList_d[tid] = minErrorPrcnt;

	}// end if (tid < nRange)

}// end function
*/

__global__ void
decode( float *pDataA, float *pDataB, float *pResult)
{

}// end function
