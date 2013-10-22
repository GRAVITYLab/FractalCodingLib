/*
 * FEL_util.h
 *
 *  Created on: Mar 14, 2012
 *      Author: abon
 */

#ifndef FEL_UTIL_H_
#define FEL_UTIL_H_

#include "ITL_util.h"
#include "ITL_vectormatrix.h"
#include "ITL_histogram.h"
//#include "ITL_distcomputer.h"
#include "ITL_entropycore.h"

template <class T>
class FEL_util
{
	enum dataTypes { SCALARDATA = 0,
					 VECTORDATA = 1,
					 TIMEVARYINGSCALARDATA = 2,
					 TIMEVARYINGVECTORDATA = 3
					};

public:

	static void
	resampleField( T* data, int dataType,
			 	   int* datasize, T** resampledData,
			 	   int* resampledsize,
			 	   float factor )
	{
		VECTOR3 f, origF, orig;
		T* temp = new T[8];
		int index = -1;

		// Calculate the size of the resampled field
		resampledsize[0] = (int)floor( datasize[0] * factor );
		resampledsize[1] = (int)floor( datasize[1] * factor );
		if( datasize[2] > 2 )
			resampledsize[2] = (int)floor( datasize[2] * factor );
		else
			resampledsize[2] = 2;

		//#ifdef DEBUG_MODE
		fprintf( stderr, "%s: %d: Field size: <%d, %d, %d>\n", __FILE__, __LINE__,
									datasize[0], datasize[1], datasize[2] );
		fprintf( stderr, "%s: %d: Downsampled field size: <%d, %d, %d>\n", __FILE__, __LINE__,
									resampledsize[0], resampledsize[1], resampledsize[2] );
		//#endif

		// Allocate memory for the resampled field
		int nVertResampledField = ITL_util<int>::prod( resampledsize, 3 );
		int nVert = datasize[0]*datasize[1]*datasize[2];
		cout << nVertResampledField << endl;
		(*resampledData) = new T[nVertResampledField];

		// Scan through the original field
		// and resample
		for( int z = 0; z<resampledsize[2]; z++ )
		{
			for( int y = 0; y<resampledsize[1]; y++ )
			{
				for( int x = 0; x<resampledsize[0]; x++ )
				{
	            	if( resampledsize[2] > 2 )
	            		origF.Set( x/factor, y/factor, z/factor );
	            	else
	            		origF.Set( x/factor, y/factor, z );

	            	orig.Set( floor( origF[0] ), floor( origF[1] ), floor( origF[2] ) );
	            	f.Set( origF[0]-orig[0], origF[1]-orig[1], origF[1]-orig[1] );

	            	//cout << x << " " << y << " " << z << " " << endl;

		            // 0
		            index = ITL_util<T>::index3DTo1D( (int)orig[0], (int)orig[1], (int)orig[2], datasize );
		            if( index < 0 || index > (nVert-1) )
		            	cout << "main index out of bound: " << index << endl;
		            //cout << index << endl;
		            temp[0] = data[index];
		            // 1
		            index = ITL_util<T>::index3DTo1D( (int)min( (int)(orig[0]+1), datasize[0]-1 ),
		            										(int)orig[1], (int)orig[2], datasize );
		            //cout << index << endl;
		            if( index < 0 || index > (nVert-1) )
		            	cout << "main index out of bound: " << index << endl;

		            temp[1] = data[index];
		            // 2
		            index = ITL_util<T>::index3DTo1D( (int)min( (int)(orig[0]+1), datasize[0]-1 ),
		            										(int)min( (int)(orig[1]+1), datasize[1]-1 ),
		            										(int)orig[2], datasize );
		            //cout << index << endl;
		            if( index < 0 || index > (nVert-1) )
		            	cout << "main index out of bound: " << index << endl;

		            temp[2] = data[index];
		            // 3
		            index = ITL_util<T>::index3DTo1D( (int)orig[0],
		            										(int)min( (int)(orig[1]+1), datasize[1]-1 ),
		            										(int)orig[2], datasize );
		            //cout << index << endl;
		            if( index < 0 || index > (nVert-1) )
		            	cout << "main index out of bound: " << index << endl;

		            temp[3] = data[index];
		            // 4
		            index = ITL_util<T>::index3DTo1D( (int)orig[0],
		            										(int) orig[1],
		            										(int)min( (int)(orig[2]+1), datasize[2]-1 ), datasize );
		            //cout << index << endl;
		            if( index < 0 || index > (nVert-1) )
		            	cout << "main index out of bound: " << index << endl;

		            temp[4] = data[index];
		            // 5
		            index = ITL_util<T>::index3DTo1D( (int)min( (int)(orig[0]+1), datasize[0]-1 ),
		            										(int)orig[1],
		            										(int)min( (int)(orig[2]+1), datasize[2]-1 ), datasize );
		            //cout << index << endl;
		            if( index < 0 || index > (nVert-1) )
		            	cout << "main index out of bound: " << index << endl;

		            temp[5] = data[index];
		            // 6
		            index = ITL_util<T>::index3DTo1D( (int)min( (int)(orig[0]+1), datasize[0]-1 ),
		            										(int)min( (int)(orig[1]+1), datasize[1]-1 ),
		            										(int)min( (int)(orig[2]+1), datasize[2]-1 ), datasize );
		            //cout << index << endl;
		            if( index < 0 || index > (nVert-1) )
		            	cout << "main index out of bound: " << index << endl;

		            temp[6] = data[index];
		            // 7
		            index = ITL_util<T>::index3DTo1D( (int)orig[0],
		            										(int)min( (int)(orig[1]+1), datasize[1]-1 ),
		            										(int)min( (int)(orig[2]+1), datasize[2]-1 ), datasize );
		            //cout << index << endl;
		            if( index < 0 || index > (nVert-1) )
		            	cout << "main index out of bound: " << index << endl;

		            temp[7] = data[index];

		            // Interpolate
		            index = ITL_util<T>::index3DTo1D( x, y, z, resampledsize );
		            if( index < 0 || index > (nVertResampledField-1) )
		            	cout << "index out of bound: " << index << endl;
		            //cout << "-> " << index << endl;
		            if( dataType == SCALARDATA || dataType == TIMEVARYINGSCALARDATA )
		            	(*resampledData)[index] = ITL_util<SCALAR>::triLinterp_scalar( (SCALAR*)temp, f );
		            else
		            	(*resampledData)[index] = ITL_util<VECTOR3>::triLinterp_vector( (VECTOR3*)temp, f );
		            //cout << "->-> " << endl;

				}
			}
		}

		delete [] temp;

	}// end function

	static void
	computeHistogramBinField( int dataType,
	  					      ITL_field_regular<T>** datafield,
	  					      ITL_field_regular<int>** binfield,
	  					      ITL_histogram* histogram, int nBin,
	  					      int* fieldSize )
	{
		assert( (*datafield) != NULL );
		assert( (*binfield) != NULL );
		SCALAR nextScalar;
		VECTOR3 nextV;
		SCALAR minV, maxV, rangeV, binWidth;

		//fprintf( stderr, "here" );
		if( dataType == SCALARDATA || dataType == TIMEVARYINGSCALARDATA )
		{
			int nV = (*datafield)->getSize();
			minV = ITL_util<T>::Min( (*datafield)->getDataFull(), nV );
			maxV = ITL_util<T>::Max( (*datafield)->getDataFull(), nV );
			rangeV = ( maxV - minV );
			binWidth  = rangeV / (float)nBin;
		}

		// Scan through each point of the histogram field
		// and convert field value to bin ID
		//fprintf( stderr, "%d %d %d\n", fieldSize[0], fieldSize[1], fieldSize[2] );
		int index1d = 0;
		int binId = 0;
		for( int z=0; z<fieldSize[2]; z++ )
		{
			for( int y=0; y<fieldSize[1]; y++ )
			{
				for( int x=0; x<fieldSize[0]; x++ )
				{
					if( dataType == SCALARDATA || dataType == TIMEVARYINGSCALARDATA )
					{
						// Get scalar at location
						nextScalar = (SCALAR)(*datafield)->getDataAt( index1d );

						// Obtain the binID corresponding to the value at this location
						binId = (int)floor( ( nextScalar - minV ) / binWidth  );
						binId = ITL_util<int>::clamp( binId, 0, nBin-1 );
						(*binfield)->setDataAt( index1d, binId );

					}
					else if( dataType == VECTORDATA || dataType == TIMEVARYINGVECTORDATA )
					{
						//fprintf( stderr, "%d %d %d\n", x ,y,z );
						// Get vector at location
						nextV = (VECTOR3)(*datafield)->getDataAt( index1d );

						// Obtain the binID corresponding to the value at this location
						binId =  histogram->get_bin_number_3D( nextV );
						binId = ITL_util<int>::clamp( binId, 0, nBin-1 );
						(*binfield)->setDataAt( index1d, binId );
					}

					// increment to the next grid vertex
					index1d += 1;

				}// end for : x
			}// end for : y
		}// end for : z

	}// end function

	static void
	computeHistogramBinField2D( int dataType,
								ITL_field_regular<T>** datafield,
								ITL_field_regular<int>** binfield,
								ITL_histogram* histogram, int nBin,
								int* fieldSize )
	{
		assert( (*datafield) != NULL );
		assert( (*binfield) != NULL );
		SCALAR nextScalar;
		VECTOR3 nextV;
		SCALAR minV, maxV, rangeV, binWidth;

		fprintf( stderr, "here 2D" );
		if( dataType == SCALARDATA || dataType == TIMEVARYINGSCALARDATA )
		{
			int nV = (*datafield)->getSize();
			minV = ITL_util<T>::Min( (*datafield)->getDataFull(), nV );
			maxV = ITL_util<T>::Max( (*datafield)->getDataFull(), nV );
			rangeV = ( maxV - minV );
			binWidth  = rangeV / (float)nBin;
		}

		// Scan through each point of the histogram field
		// and convert field value to bin ID
		int index1d = 0;
		int binId = 0;
		for( int z=0; z<fieldSize[2]; z++ )
		{
			for( int y=0; y<fieldSize[1]; y++ )
			{
				for( int x=0; x<fieldSize[0]; x++ )
				{
					if( dataType == SCALARDATA || dataType == TIMEVARYINGSCALARDATA )
					{
						// Get scalar at location
						nextScalar = (SCALAR)(*datafield)->getDataAt( index1d );

						// Obtain the binID corresponding to the value at this location
						binId = (int)floor( ( nextScalar - minV ) / binWidth  );
						binId = ITL_util<int>::clamp( binId, 0, nBin-1 );
						(*binfield)->setDataAt( index1d, binId );

					}
					else if( dataType == VECTORDATA || dataType == TIMEVARYINGVECTORDATA )
					{
						// Get vector at location
						nextV = (VECTOR3)(*datafield)->getDataAt( index1d );

						// Obtain the binID corresponding to the value at this location
						binId = histogram->get_bin_number_2D( nextV, nBin );
						binId = ITL_util<int>::clamp( binId, 0, nBin-1 );
						(*binfield)->setDataAt( index1d, binId );

					}

					// increment to the next grid vertex
					index1d += 1;

				}// end for : x
			}// end for : y
		}// end for : z

	}// end function

	static void
	computeHistogramBinField( int dataType,
	  					      ITL_field_regular<T>*** datafield,
	  					      ITL_field_regular<int>** binfield,
	  					      ITL_histogram* histogram, int nBin,
	  					      int* fieldSize )
	{
		assert( (*datafield) != NULL );
		assert( (*binfield) != NULL );
		SCALAR nextScalar;
		VECTOR3 nextV;
		SCALAR minV, maxV, rangeV, binWidth;
		ITL_field_regular<T>** tmp = (*datafield);

		if( dataType == SCALARDATA || dataType == TIMEVARYINGSCALARDATA )
		{
			int nV = (*tmp)->getSize();
			minV = ITL_util<T>::Min( (*tmp)->getDataFull(), nV );
			maxV = ITL_util<T>::Max( (*tmp)->getDataFull(), nV );
			rangeV = ( maxV - minV );
			binWidth  = rangeV / (float)nBin;
		}

		// Scan through each point of the histogram field
		// and convert field value to bin ID
		int index1d = 0;
		int binId = 0;
		for( int z=0; z<fieldSize[2]; z++ )
		{
			for( int y=0; y<fieldSize[1]; y++ )
			{
				for( int x=0; x<fieldSize[0]; x++ )
				{
					if( dataType == SCALARDATA || dataType == TIMEVARYINGSCALARDATA )
					{
						// Get scalar at location
						nextScalar = (SCALAR)(*tmp)->getDataAt( index1d );

						// Obtain the binID corresponding to the value at this location
						binId = (int)floor( ( nextScalar - minV ) / binWidth  );
						binId = ITL_util<int>::clamp( binId, 0, nBin-1 );
						(*binfield)->setDataAt( index1d, binId );

					}
					else if( dataType == VECTORDATA || dataType == TIMEVARYINGVECTORDATA )
					{
						// Get vector at location
						nextV = (VECTOR3)(*tmp)->getDataAt( index1d );

						// Obtain the binID corresponding to the value at this location
						binId =  histogram->get_bin_number_3D( nextV );
						binId = ITL_util<int>::clamp( binId, 0, nBin-1 );
						(*binfield)->setDataAt( index1d, binId );
					}

					// increment to the next grid vertex
					index1d += 1;

				}// end for : x
			}// end for : y
		}// end for : z

	}// end function

	static void
	computeHistogramBinField2D( int dataType,
								ITL_field_regular<T>*** datafield,
								ITL_field_regular<int>** binfield,
								ITL_histogram* histogram, int nBin,
								int* fieldSize )
	{
		assert( (*datafield) != NULL );
		assert( (*binfield) != NULL );
		SCALAR nextScalar;
		VECTOR3 nextV;
		SCALAR minV, maxV, rangeV, binWidth;
		ITL_field_regular<T>** tmp = (*datafield);

		if( dataType == SCALARDATA || dataType == TIMEVARYINGSCALARDATA )
		{
			int nV = (*tmp)->getSize();
			minV = ITL_util<T>::Min( (*tmp)->getDataFull(), nV );
			maxV = ITL_util<T>::Max( (*tmp)->getDataFull(), nV );
			rangeV = ( maxV - minV );
			binWidth  = rangeV / (float)nBin;
		}

		// Scan through each point of the histogram field
		// and convert field value to bin ID
		int index1d = 0;
		int binId = 0;
		for( int z=0; z<fieldSize[2]; z++ )
		{
			for( int y=0; y<fieldSize[1]; y++ )
			{
				for( int x=0; x<fieldSize[0]; x++ )
				{
					if( dataType == SCALARDATA || dataType == TIMEVARYINGSCALARDATA )
					{
						// Get scalar at location
						nextScalar = (SCALAR)(*tmp)->getDataAt( index1d );

						// Obtain the binID corresponding to the value at this location
						binId = (int)floor( ( nextScalar - minV ) / binWidth  );
						binId = ITL_util<int>::clamp( binId, 0, nBin-1 );
						(*binfield)->setDataAt( index1d, binId );

					}
					else if( dataType == VECTORDATA || dataType == TIMEVARYINGVECTORDATA )
					{
						// Get vector at location
						nextV = (VECTOR3)(*tmp)->getDataAt( index1d );

						// Obtain the binID corresponding to the value at this location
						binId = histogram->get_bin_number_2D( nextV, nBin );
						binId = ITL_util<int>::clamp( binId, 0, nBin-1 );
						(*binfield)->setDataAt( index1d, binId );

					}

					// increment to the next grid vertex
					index1d += 1;

				}// end for : x
			}// end for : y
		}// end for : z

	}// end function

	static void
	computeHistogramFromBlock( T *blockData, int nPoint, ITL_histogram* histogram, int nBin, double* freqList )
	{
		memset( freqList, 0, sizeof(double)*nBin );

		for( int iP=0; iP<nPoint; iP++ )
		{
			VECTOR3 nextV = blockData[iP];
			int binId = histogram->get_bin_number_2D( nextV, nBin );
			freqList[binId] ++;
		}
		for( int i=0; i<nBin; i++ )
			freqList[i] = freqList[i] / (double)nPoint;

	}// end function


	static void
	computeHistogramFromBinField( int *binData, int nPoint, int nBin, double* freqList )
	{
		assert( binData != NULL );
		assert( freqList != NULL );

		memset( freqList, 0.0, sizeof(double)*nBin );

		int binID = -1;
		for( int iP=0; iP<nPoint; iP++ )
		{
			binID = binData[iP];

			freqList[binID] ++;
		}

		for( int i=0; i<nBin; i++ )
			freqList[i] = freqList[i] / (double)nPoint;

	}// End function

	static void
	computeHistogramFrequencies( float* data, float* dataRange, double* freqList, int nPoint, int nbin )
	{
		int binId = 0;
		float range = ( dataRange[1] - dataRange[0] );
		float binWidth = range / (float)nbin;

		memset( freqList, 0.0, sizeof(double)*nbin );

		// Scan through data
		for( int i=0; i<nPoint; i++ )
		{
			binId = (int)floor( (data[i] - dataRange[0] ) / binWidth );
			if( binId > nbin-1 )
				binId = nbin-1;

			freqList[ binId ] ++;
		}

		for( int i=0; i<nbin; i++ )
			freqList[i] /= (float)nPoint;


	}// end function


	static float
	computeEntropy_TimeHistogram( double* tvFreqList, int nBin, int nTimeStep, bool isNormalize )
	{
		float sum = 0;
		for( int t=0; t<nTimeStep; t++ )
		{
			sum += (float)ITL_entropycore::computeEntropy_HistogramBased2( tvFreqList + t*nBin, nBin, false );
		}

		return ( sum / (float)nTimeStep );

	}// end function

	/**
	 * Apply various transformations to the domain distribution
	 * and compare to the range distribution
	 */
	static double
	matchRangeToDomainByDistribution( int distanceMeasureType,
									  double* rangeFreqList,
									  double* domainFreqList,
									  int nBin,
						 	 	 	  int* rotateAmount,
						 	 	 	  bool* isReflected )
	{
		double* domainFreqListCopy = new double[nBin*2];
		double* domainFreqListRev = new double[nBin*2];
		double* buffer2 = new double[nBin-1];
		double matchErrorPrcnt = 100;
		double minErrorPrcnt = 100;
		int optRot;
		bool optRef;
		double buffer;
		int nRotLimit = 10;

		#ifdef DEBUG_MODE
		double m = ITL_util<double>::Min( rangeFreqList, nBin );
		double M = ITL_util<double>::Max( rangeFreqList, nBin );
		fprintf( stderr, "Range distribution range: %g, %g\n", m, M );
		#endif

		// Circular shift
		// (nBin-1) circular shifts possible
		memcpy( domainFreqListCopy, domainFreqList, sizeof(double)*nBin );
		memcpy( domainFreqListCopy+nBin, domainFreqList, sizeof(double)*nBin );

		//for( int iT = 0; iT<nBin-1; iT++ )
		for( int iT = 0; iT<nRotLimit-1; iT++ )
		{
			// Compare with range distribution
			matchErrorPrcnt = FEL_util<double>:: computeHistogramMatchError( distanceMeasureType,
																			 rangeFreqList, domainFreqListCopy+iT,
																			 nBin );
			if( iT == 0 )
					minErrorPrcnt = matchErrorPrcnt;

			//fprintf( stderr, "%d: Match/min error: %g %g\n", iT, matchErrorPrcnt, minErrorPrcnt );
			if( matchErrorPrcnt < minErrorPrcnt )
			{
				//fprintf( stderr, "updating\n" );
				minErrorPrcnt = matchErrorPrcnt;
				//*optTransformationMatrix = trManager.rotMatrixArray[iT];
				optRot = nBin - iT;
				optRef = false;
			}
		}

		// Reflect distribution
		for( int i=0; i<nBin; i++ )
			domainFreqListRev[i] = domainFreqList[nBin-i-1];
		memcpy( domainFreqListRev+nBin, domainFreqListRev, sizeof(double)*nBin );

		//for( int iT = 0; iT<nBin-1; iT++ )
		for( int iT = 0; iT<nRotLimit-1; iT++ )
		{
			// Compare with range distribution
			matchErrorPrcnt = FEL_util<double>:: computeHistogramMatchError( distanceMeasureType,
																			 rangeFreqList, domainFreqListRev+iT,
																			 nBin );

			//fprintf( stderr, "%d: Match/min error: %g %g\n", iT, matchErrorPrcnt, minErrorPrcnt );
			if( matchErrorPrcnt < minErrorPrcnt )
			{
				minErrorPrcnt = matchErrorPrcnt;
				optRot = nBin - iT;
				optRef = false;
			}
		}

		(*rotateAmount) = optRot;
		(*isReflected) = optRef;
		//fprintf( stderr, "Optimal transformation: %d, %d\n", optRot, (int)optRef );

		// Free temporary resources
		delete [] domainFreqListCopy;
		delete [] domainFreqListRev;
		delete [] buffer2;

		return minErrorPrcnt;

	}// End function

	/**
	 * Apply various transformations to the domain distribution
	 * and compare to the range distribution
	 */
	static double
	matchRangeToDomainByDistribution_Fast( int distanceMeasureType,
									  	  	  	double* rangeFreqList,
									  	  	  	long domainIndex,
									  	  	  	double** dBuffer, double** dBufferRev,
									  	  	  	int nBin, int* rotateAmount, bool* isReflected )
	{
		double* buffer2 = new double[nBin-1];
		double matchErrorPrcnt = 100;
		double minErrorPrcnt = 100;
		int optRot;
		bool optRef;
		double buffer;
		int nRotLimit = floor( 0.25*nBin);

		#ifdef DEBUG_MODE
		double m = ITL_util<double>::Min( rangeFreqList, nBin );
		double M = ITL_util<double>::Max( rangeFreqList, nBin );
		fprintf( stderr, "Range distribution range: %g, %g\n", m, M );
		#endif

		//for( int iT = 0; iT<nBin-1; iT++ )
		for( int iT = 0; iT<nRotLimit-1; iT++ )
		{
			// Compare with range distribution
			matchErrorPrcnt = FEL_util<double>:: computeHistogramMatchError( distanceMeasureType,
																			 rangeFreqList,
																			 *(dBuffer + domainIndex) + iT,
																			 nBin );
			if( iT == 0 )
					minErrorPrcnt = matchErrorPrcnt;

			//fprintf( stderr, "%d: Match/min error: %g %g\n", iT, matchErrorPrcnt, minErrorPrcnt );
			if( matchErrorPrcnt < minErrorPrcnt )
			{
				//fprintf( stderr, "updating\n" );
				minErrorPrcnt = matchErrorPrcnt;
				//*optTransformationMatrix = trManager.rotMatrixArray[iT];
				optRot = nBin - iT;
				optRef = false;
			}
		}

		//for( int iT = 0; iT<nBin-1; iT++ )
		for( int iT = 0; iT<nRotLimit-1; iT++ )
		{
			// Compare with range distribution
			matchErrorPrcnt = FEL_util<double>:: computeHistogramMatchError( distanceMeasureType,
																			 rangeFreqList,
																			 *(dBufferRev + domainIndex) + iT,
																			 nBin );

			//fprintf( stderr, "%d: Match/min error: %g %g\n", iT, matchErrorPrcnt, minErrorPrcnt );
			if( matchErrorPrcnt < minErrorPrcnt )
			{
				minErrorPrcnt = matchErrorPrcnt;
				optRot = nBin - iT;
				optRef = false;
			}
		}

		(*rotateAmount) = optRot;
		(*isReflected) = optRef;
		//fprintf( stderr, "Optimal transformation: %d, %d\n", optRot, (int)optRef );

		// Free temporary resources
		delete [] buffer2;

		return minErrorPrcnt;

	}// End function

	static void
	compareDistributionWithTemplate( T* iHist,
									 	 T* templateHist,
									 	 int nBin,
									 	 T* error,
									 	 T* errorList,
									 	 int* r_id,
									 	 bool* ref_id )
	{
		T matchError = 0;
		T transformedTemplate[nBin];
		T locallyTransformedTemplate[nBin];

		//memcpy( transformedTemplate, templateHist, sizeof(T)*nBin );

		// Find peak of distribution
		int iPeak = findPeakLocation( iHist, nBin );
		//printHistogram( iHist, nBin );
		//fprintf( stderr, "My peak is at %d\n", iPeak );

		// Find peak of template
		//printHistogram( templateHist, nBin );
		int tPeak = findPeakLocation( templateHist, nBin );
		//fprintf( stderr, "Template's peak is at %d\n", tPeak );

		// Transform template
		// so that after transformation,
		// its peak matches the IH's peak
		//shiftDistribution_zerofill( templateHist, transformedTemplate, nBin, (iPeak-tPeak) );
		shiftDistribution_rotate( templateHist, transformedTemplate, nBin, (iPeak-tPeak) );
		//printHistogram( transformedTemplate, nBin );

		for( int i=-2; i<=2; i++ )
		{
			//shiftDistribution_zerofill( transformedTemplate, locallyTransformedTemplate, nBin, i );
			shiftDistribution_rotate( transformedTemplate, locallyTransformedTemplate, nBin, i );

			// Compute error
			matchError = getL2Norm( iHist, locallyTransformedTemplate, nBin );
			//fprintf( stderr, "i=%d, Match error: %d ...\n", i, matchError );
			//printHistogram( locallyTransformedTemplate, nBin );

			// Update min error
			if( i == -2 )
			{
				(*error) = matchError;
				(*r_id) = (iPeak-tPeak) + i;
				(*ref_id) = false;

				for( int iB=0; iB<nBin; iB++ )
					errorList[iB] = ( iHist[iB] - locallyTransformedTemplate[iB] );
			}
			if( matchError < (*error) )
			{
				(*error) = matchError;
				(*r_id) = (iPeak-tPeak) + i; // Compound transformation
				(*ref_id) = false;

				for( int iB=0; iB<nBin; iB++ )
					errorList[iB] = ( iHist[iB] - locallyTransformedTemplate[iB] );

			}
		}// end for

	}// end function

	static void
	transformDistribution( double* freqList, int nBin, int nRot, bool isRef )
	{
		double buffer[nBin-nRot];
		double buffer2[nRot];

		if( nRot == 0 && isRef == false )
		{
			//fprintf( stderr, "trivial rot/ref\n" );
			return;
		}

		// Create buffer
		double freqListTemp[nBin];

		// Reflect, if needed
		if( isRef )
		{
			for( int i=0; i<nBin; i++ )
				freqListTemp[i] = freqList[nBin-i-1];
		}
		else
			memcpy( freqListTemp, freqList, sizeof(double)*nBin );

		if( nRot > 0 )
		{
			memcpy( buffer, freqListTemp, sizeof(double)*(nBin-nRot) );
			memcpy( buffer2, freqListTemp+(nBin-nRot), sizeof(double)*nRot );
			memcpy( freqListTemp, buffer2, sizeof(double)*nRot );
			memcpy( freqListTemp+nRot, buffer, sizeof(double)*(nBin-nRot) );
		}

		memcpy( freqList, freqListTemp, sizeof(double)*nBin );

	}// End function

	static void
	rotateDistribution( T* freqList, int nBin, int nShift )
	{
		T freqListTemp[nBin];
		T buffer[nBin-nShift];
		T buffer2[nBin];

		if( nShift == 0 )
			return;

		memcpy( freqListTemp, freqList, sizeof(T)*nBin );

		if( nShift > 0 )
		{
			memcpy( buffer, freqListTemp, sizeof(T)*(nBin-nShift) );
			memcpy( buffer2, freqListTemp+(nBin-nShift), sizeof(T)*nShift );
			memcpy( freqListTemp, buffer2, sizeof(T)*nShift );
			memcpy( freqListTemp+nShift, buffer, sizeof(T)*(nBin-nShift) );
		}

		memcpy( freqList, freqListTemp, sizeof(T)*nBin );

	}// End function

	static void
	shiftDistribution_zerofill( T* freqListSrc, T* freqListDst,
									int nBin, int nShift )
	{
		T freqListTemp[nBin];
		T buffer[nBin-nShift];

		if( nShift == 0 )
		{
			memcpy( freqListDst, freqListSrc, sizeof(T)*nBin );
			return;
		}

		memcpy( freqListTemp, freqListSrc, sizeof(long)*nBin );

		// Left - clockwise
		if( nShift < 0 )
		{
			memcpy( freqListTemp, freqListSrc+abs(nShift), sizeof(T)*( nBin-abs(nShift) ) );
			memset( freqListTemp + (nBin-abs(nShift) ), 0, sizeof(T)*abs(nShift) );
		}
		// Right - anti-clockwise
		else if( nShift > 0 )
		{
			memcpy( freqListTemp + nShift, freqListSrc, sizeof(T)*(nBin-nShift) );
			memset( freqListTemp, 0, sizeof(T)*nShift );
		}

		memcpy( freqListDst, freqListTemp, sizeof(T)*nBin );

	}// End function

	static void
	shiftDistribution_rotate( T* freqListSrc, T* freqListDst,
									int nBin, int nShift )
	{
		T freqListTemp[nBin];
		T buffer[nBin-nShift];

		if( nShift == 0 )
		{
			memcpy( freqListDst, freqListSrc, sizeof(T)*nBin );
			return;
		}

		memcpy( freqListTemp, freqListSrc, sizeof(long)*nBin );

		// Left - clockwise
		if( nShift < 0 )
		{
			memcpy( freqListTemp, freqListSrc+abs(nShift), sizeof(T)*( nBin-abs(nShift) ) );
			int start = nBin-abs(nShift);
			for( int i=0; i<abs(nShift); i++ )
				freqListTemp[start+i] = freqListSrc[abs(nShift)-1-i];
		}
		// Right - anti-clockwise
		else if( nShift > 0 )
		{
			memcpy( freqListTemp + nShift, freqListSrc, sizeof(T)*(nBin-nShift) );
			for( int i=0; i<nShift; i++ )
				freqListTemp[i] = freqListSrc[nBin-1-i];
		}

		memcpy( freqListDst, freqListTemp, sizeof(T)*nBin );

	}// End function


	static double
	computeHistogramMean( double* freqList, int nBin )
	{
		double mu = 0;
		double binwidth = 1 / (double)nBin;

		for( int i=0; i<nBin; i++ )
		{
			mu += ( freqList[i] * ((i+0.5)*binwidth) );
		}

		return mu;
	}

	static double
	computeHistogramMean( double* bincenterList, double* freqList, int nBin )
	{
		double mu = 0.0, w = 0.0;

		for( int i=0; i<nBin; i++ )
		{
			mu += ( freqList[i] * bincenterList[i] );
			w += bincenterList[i];
		}

		return ( mu / w );
	}

	static double
	computeHistogramSD( double* freqList, int nBin )
	{
		double mu = FEL_util<double>::computeHistogramMean( freqList, nBin );
		double binwidth = 1 / (double)nBin;
		double sumsq = 0;

		for( int i=0; i<nBin; i++ )
		{
			sumsq += freqList[i] * ( (i+0.5)*binwidth - mu ) * ( (i+0.5)*binwidth - mu );
		}

		return sqrt( sumsq );
	}

	static double
	computeHistogramMeanError( double* f1List, double* f2List, int nBin )
	{
		return abs( FEL_util<double>::computeHistogramMean( f1List, nBin ) -
					FEL_util<double>::computeHistogramMean( f2List, nBin ) );
	}

	static double
	computeHistogramSDError( double* f1List, double* f2List, int nBin )
	{
		return abs( FEL_util<double>::computeHistogramSD( f1List, nBin ) -
					FEL_util<double>::computeHistogramSD( f2List, nBin ) );
	}

	static double
	computeHistogramMatchError( int distanceMeasureType,
								double* rangeDataHist, double* trDomainDataHist,
								int nBin )
	{
		double dist = 0;
		//fprintf( stderr, "%d", distanceMeasureType );
		if( distanceMeasureType == 0 ) 			dist = FEL_util<double>::computeL1( rangeDataHist, trDomainDataHist, nBin );
		else if( distanceMeasureType == 1 ) 	dist = FEL_util<double>::computeL2( rangeDataHist, trDomainDataHist, nBin );
		else if( distanceMeasureType == 2 ) 	dist = FEL_util<double>::computeHistIntersection( rangeDataHist, trDomainDataHist, nBin );
		else if( distanceMeasureType == 3 ) 	dist = FEL_util<double>::computeEMD1D( rangeDataHist, trDomainDataHist, nBin );
		else if( distanceMeasureType == 4 ) 	dist = FEL_util<double>::computeChiSquare( rangeDataHist, trDomainDataHist, nBin );
		//else if( distMeasureType == 4 ) 	dist = ITL_distcomputer::computeKLD( rangeDataHist, trDomainDataHist, nBin );
		//else if( distMeasureType == 5 ) 	dist = ITL_distcomputer::computeJSD( rangeDataHist, trDomainDataHist, nBin );
		//fprintf( stderr, "%g", dist );

		return dist;

	}// End function

	static double
	computeL1( double *dist1, double *dist2, int nBin )
	{
		assert( dist1 != NULL );
		assert( dist2 != NULL );

		double dk = 0;
		double l1 = 0;
		for(int iB=0; iB<nBin; iB++ )
		{
			//printf( "%d: %f %f\n", iB, dist1[iB], dist2[iB] );
			dk = dist1[iB]-dist2[iB];
			l1 += abs(dk);
		}

		return l1;
	}

	static double
	computeL2( double *dist1, double *dist2, int nBin )
	{
		assert( dist1 != NULL );
		assert( dist2 != NULL );

		double dk = 0;
		double l2 = 0;
		for(int iB=0; iB<nBin; iB++ )
		{
			//printf( "%d: %f %f\n", iB, dist1[iB], dist2[iB] );
			dk = dist1[iB]-dist2[iB];
			l2 += (dk*dk);
		}

		return sqrt(l2);
	}

	static double
	computeHistIntersection( double *dist1, double *dist2, int nBin )
	{
		assert( dist1 != NULL );
		assert( dist2 != NULL );

		double l = 0;
		for( int iB=0; iB<nBin; iB++ )
		{
			l += std::min( dist1[iB], dist2[iB] );
		}

		return (1-l);
	}

	static double
	computeChiSquare( double *dist1, double *dist2, int nBin )
	{
		assert( dist1 != NULL );
		assert( dist2 != NULL );

		double l = 0;
		for( int iB=0; iB<nBin; iB++ )
		{
			l += ( ( dist1[iB] - dist2[iB] )*( dist1[iB] - dist2[iB] ) ) / ( ( dist1[iB] + dist2[iB] )/2.0f ) ;
		}

		return l;
	}

	static double
	computeKLD( double *dist1, double *dist2, int nBin )
	{
		return 0.0f;
	}

	static double
	doublecomputeJSD( double *dist1, double *dist2, int nBin )
	{
		return 0.0f;
	}

	static double
	computeEMD1D( double *dist1, double *dist2, int nBin )
	{
		assert( dist1 != NULL );
		assert( dist2 != NULL );
		#if defined( _WIN32 ) || defined( _WIN64 )
			float* cDist1 = new float[nBin];
			float* cDist2 = new float[nBin];
		#else
			double cDist1[nBin];
			double cDist2[nBin];
		#endif

		double dk = 0;
		double l3 = 0;
		cDist1[0] = cDist2[0] = 0;
		for(int iB=0; iB<nBin; iB++ )
		{
			//printf( "%d: %f %f\n", iB, dist1[iB], dist2[iB] );

			if( iB == 0 )
			{
				cDist1[iB] = dist1[iB];
				cDist2[iB] = dist2[iB];
			}
			else
			{
				cDist1[iB] = cDist1[iB-1] + dist1[iB];
				cDist2[iB] = cDist2[iB-1] + dist2[iB];
			}

			dk = cDist1[iB]-cDist2[iB];
			l3 += abs(dk);
			//printf( "l3: %g\n", l3 );
		}

		#if defined( _WIN32 ) || defined( _WIN64 )
			delete [] cDist1;
			delete [] cDist2;
		#endif

		return l3;
	}// End function


	static double
	computeTimeHistogramMatchError( int distanceMeasureType,
									double* rangeDataTvHist, double* trDomainDataTvHist,
									int nBin, int nTime )
	{
		double dist;
		if( distanceMeasureType == 0 ) 			dist = FEL_util<T>::computeTimeHistL1( rangeDataTvHist, trDomainDataTvHist, nBin, nTime );
		else if( distanceMeasureType == 1 ) 	dist = FEL_util<T>::computeTimeHistL2( rangeDataTvHist, trDomainDataTvHist, nBin, nTime );
		else if( distanceMeasureType == 2 ) 	dist = FEL_util<T>::computeTimeHistIntersection( rangeDataTvHist, trDomainDataTvHist, nBin, nTime );
		//else if( distanceMeasureType == 3 ) 	dist = ITL_distcomputer::computeEMD1D( rangeDataHist, trDomainDataHist, nBin );
		//else if( distMeasureType == 4 ) 	dist = ITL_distcomputer::computeKLD( rangeDataHist, trDomainDataHist, nBin );
		//else if( distMeasureType == 5 ) 	dist = ITL_distcomputer::computeJSD( rangeDataHist, trDomainDataHist, nBin );

		return dist;

	}// End function

	static double
	computeTimeHistL1( double* rangeDataTvHist, double* trDomainDataTvHist,
					  int nBin, int nTime )
	{
		fprintf( stderr, "Not implemented yet\n" );
		return 0.0;
	}

	static double
	computeTimeHistL2( double* rangeDataTvHist, double* trDomainDataTvHist,
			int nBin, int nTime )
	{
		fprintf( stderr, "Not implemented yet\n" );
		return 0.0;
	}

	static double
	computeTimeHistIntersection( double* rangeDataTvHist, double* trDomainDataTvHist,
								 int nBin, int nTime )
	{
		double distSum = 0;
		for( int t=0; t<nTime; t++ )
		{
			distSum += FEL_util<double>::computeHistIntersection( rangeDataTvHist + t*nBin,
																  trDomainDataTvHist + t*nBin,
																  nBin );
		}

		return distSum / (float)nTime;
	}


	static bool
	isOverlapping( T* a, T* b )
	{
		#ifdef DEBUG_MODE
		cout << a[0] << " " << a[1] << " " << a[2] << " " << a[3] << " " << a[4] << " " << a[5] <<  endl;
		cout << b[0] << " " << b[1] << " " << b[2] << " " << b[3] << " " << b[4] << " " << b[5] << endl;
		#endif

		if( a[3] <= b[0] || b[3] <= a[0] ||
			a[4] <= b[1] || b[4] <= a[1] ||
			a[5] <= b[2] || b[5] <= a[2] )
		{
			//cout << "overlap: no" << endl;
			return false;
		}

		//cout << "overlap: yes" << endl;
		return true;
	};

	static int
	getOverlapSize( T* a, T* b )
	{
		#ifdef DEBUG_MODE
		cout << a[0] << " " << a[1] << " " << a[2] << " " << a[3] << " " << a[4] << " " << a[5] <<  endl;
		cout << b[0] << " " << b[1] << " " << b[2] << " " << b[3] << " " << b[4] << " " << b[5] << endl;
		#endif

		int n = 0;

		if( a[3] <= b[0] || b[3] <= a[0] ||
			a[4] <= b[1] || b[4] <= a[1] ||
			a[5] <= b[2] || b[5] <= a[2] )
		{
			//fprintf( stderr, "overlap: n" );
			n = 0;
		}
		else
		{
			//fprintf( stderr, "overlap: y" );
			n = max( a[3]-b[0]+1, b[3]-a[0] + 1 ) * max( a[4]-b[1]+1, b[4]-a[1] + 1 ) * max( a[5]-b[2]+1, b[5]-a[2] + 1 );
		}

		return n;
	};

	static void
	scaleDistribution( double* fList, int nBin, double scale )
	{
		double sum = 0;
		for( int i=0; i<nBin; i++ )
		{
			fList[i] = fList[i]*scale;
			sum += fList[i];
		}
		//for( int i=0; i<nBin; i++ )
		//{
		//	fList[i] = fList[i] / sum;

		//}

	}

	static double
	getL2Norm( T* a, T* b, int len )
	{
		T sqSum = 0;
		for( int i=0; i<len; i++ )
			sqSum += ( a[i] - b[i] ) * ( a[i] - b[i] );

		return sqrt( sqSum / len );
	}

	static int
	findPeakLocation( T* iHist, int nbin )
	{
		int max_id = 0;
		T maxFreq = 0;

		for( int i=0; i<nbin; i++ )
		{
			if( iHist[i] > maxFreq )
			{
				max_id = i;
				maxFreq = iHist[i];
			}
		}

		return max_id;

	}// end function

};

#endif
/* FEL_UTIL_H_ */
