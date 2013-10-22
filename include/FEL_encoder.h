/*
 * FEL_encoder.h
 *
 *  Created on: Feb 15, 2012
 *      Author: abon
 */

#ifndef FEL_ENCODER_H_
#define FEL_ENCODER_H_

#include <list>

#include "ITL_distcomputer.h"
#include "ITL_distribution.h"
#include "ITL_histogram.h"
#include "ITL_entropycore.h"
#include "ITL_field_regular.h"
#include "ITL_partition_quadtree.h"
#include "ITL_partition_octree.h"

#include "FEL_util.h"
#include "FEL_range.h"
#include "FEL_domain.h"
#include "FEL_transformation.h"
#include "FEL_encodedhistogram.h"
#include "FEL_spacetree.h"
#include "FEL_segmenttree.h"
#include "FEL_matlab.h"

template <class T>
class FEL_encoder
{
	enum dataTypes { SCALARDATA = 0,
					 VECTORDATA = 1 };

	enum encodingTypes { REGULAR = 0,
						 QUADTREE = 1,
						 OCTREE = 2 };

	enum transformTypes { ROTATION = 0,
						  REFLECTION = 1 };

	enum domainCreationTypes{ SINGLERES = 0,
							  MULTIRES = 1
							};

	int dataType;

	int nDim;
	int nTotalRange, nLeafRange;
	int nTotalDomain;
	int nUnMatchedRange;
	int nTotalHighErrorBin;


	float minDomainEntropy, maxDomainEntropy;

	bool isSharingNeighbor;

	ITL_histogram *histogram;

	list< FEL_range<T> > rangeBlockList;
	list< FEL_domain<T> > domainList;

	int globalRangeIndex;

	FEL_spacetree<T> rangeBlockTree;
	ITL_partition_quadtree<T> quadtreePartitioner;
	ITL_partition_octree<T> octreePartitioner;

	FEL_segmenttree<T> domainTree;

	FEL_transformation trManager;

public:

	list<ITL_distribution> salientRangeDistributionList;

	bool isDomainPruningOn;

	//double resultRegionLimitList[175]; 	// 25 * 7
	float resultRegionLimitList[175]; 		// 25 * 7

	// Histogram properties
	int histMappingType;
	int nBin;

	// Domain creation parameters
	int domainCreationType;

	// Encoding parameters
	int encodingType;
	int distanceMeasureType;
	int smallestDomainSize[4];
	int smallestRangeSize[4];
	float matchErrorThreshold;
	float entropyThreshold;

	// Performance
	double domainCreationTime, domainOrgTime, encodeTime, searchTime;

	// Accuracy Analysis
	bool matlab_plotting_on;

	// Plotting
	FEL_matlab* plotter;

	ITL_field_regular<T>** inputDataField;
	ITL_field_regular<int>* inputBinField;

	ITL_field_regular<T>* downsampledDataField;
	ITL_field_regular<int>* downsampledBinField;
	ITL_field_regular<T>* downsampledDataFieldArray;
	ITL_field_regular<int>* downsampledBinFieldArray;

	// Constructor
	FEL_encoder()
	{
		// Set default values to some parameters
		dataType = SCALARDATA;
		encodingType = REGULAR;
		histMappingType = 2;
		domainCreationType = MULTIRES;
		distanceMeasureType = 0;
		isDomainPruningOn = false;
		//matlab_plotting_and_analysis_on = false;

		// Initialize pointers to NULL
		inputDataField = NULL;
		inputBinField = NULL;
		histogram = NULL;
		//plotter = NULL;

	}// end constructor

	FEL_encoder( int datatype,
				 ITL_field_regular<T>** infield,
				 int *domainsize, int *rangesize )
	{
		// Get data type
		dataType = datatype;

		// Get data
		inputDataField = infield;
		nDim = (*inputDataField)->getNumDim();

		// Get data size
		float lowF[nDim];
		float highF[nDim];
		int dim[nDim];
		//inputDataField->getBounds( lowF, highF );
		//inputDataField->getSize( dim );
		(*inputDataField)->getBounds( lowF, highF );
		(*inputDataField)->getSize( dim );

		// Create histogram bin field
		inputBinField = new ITL_field_regular<int>();
		inputBinField->initialize( nDim, lowF, highF );
		if( histMappingType == 2 )
			FEL_util<T>::computeHistogramBinField2D( dataType,
													 &inputDataField, &inputBinField,
													 histogram, nBin,
													 dim );
		else
			FEL_util<T>::computeHistogramBinField( dataType,
												   &inputDataField, &inputBinField,
												   histogram, nBin,
												   dim );

		// Record smallest domain and range sizes
		for( int i=0; i<nDim; i++ )
		{
			smallestDomainSize[i] = domainsize[i];
			smallestRangeSize[i] = rangesize[i];
		}

		// Initialize some other parameters
		nTotalDomain = 0;
		globalRangeIndex = 0;
		isSharingNeighbor = false;
		isDomainPruningOn = false;
		//matlab_plotting_and_analysis_on = false;

		// Initialize pointers to NULL
		histogram = NULL;
		//plotter = NULL;

	}// End constructor

	~FEL_encoder()
	{
		domainList.clear();
		rangeBlockList.clear();

		if( inputBinField != NULL )	delete inputBinField;
		if( downsampledDataField != NULL ) delete downsampledDataField;
		if( downsampledBinField != NULL ) delete downsampledBinField;
		if( downsampledDataFieldArray != NULL ) delete [] downsampledDataFieldArray;
		if( downsampledBinFieldArray != NULL ) delete [] downsampledBinFieldArray;

	}// End destructor

	void
	initialize( int datatype,
				ITL_field_regular<T>** infield,
	  	 	    int *domainsize, int *rangesize )
	{
		// Get data type
		dataType = datatype;

		// Get data
		inputDataField = infield;
		//nDim = inputDataField->getNumDim();
		nDim = (*inputDataField)->getNumDim();
		cout << "1" << endl;

		// Get data size
		float lowF[nDim];
		float highF[nDim];
		int dim[nDim];
		//inputDataField->getBounds( lowF, highF );
		//inputDataField->getSize( dim );
		(*inputDataField)->getBounds( lowF, highF );
		(*inputDataField)->getSize( dim );
		cout << "2" << endl;

		// Create histogram bin field
		inputBinField = new ITL_field_regular<int>();
		cout << "3" << endl;

		inputBinField->initialize( nDim, lowF, highF );
		cout << "4" << endl;

		if( histMappingType == 2 )
		{
			cout << "2d" << endl;
			FEL_util<T>::computeHistogramBinField2D( dataType,
													 &inputDataField, &inputBinField,
													 histogram, nBin,
													 dim );
		}
		else
		{
			cout << "3d" << nBin << endl;
			assert( histogram != NULL );
			assert( inputDataField != NULL );
			assert( inputBinField != NULL );
			FEL_util<T>::computeHistogramBinField( dataType,
												   &inputDataField, &inputBinField,
												   histogram, nBin,
												   dim );
		}

		cout << "5" << endl;

		// Record smallest domain and range sizes
		for( int i=0; i<nDim; i++ )
		{
			smallestDomainSize[i] = domainsize[i];
			smallestRangeSize[i] = rangesize[i];
		}

		// Initialize transformation class
		//trManager.initTransformationMatrices( 16 );

		// Initialize some other variables
		nTotalDomain = 0;
		globalRangeIndex = 0;
		isSharingNeighbor = false;
		//matlab_plotting_and_analysis_on = false;
		//plotter = NULL;

	}// End constructor

	void createDomains()
	{
		assert( inputDataField != NULL );
		assert( inputBinField != NULL );

		clock_t starttime = ITL_util<float>::startTimer();

		if( domainCreationType == SINGLERES )
			createDomainsMultiSize();
		else if( domainCreationType == MULTIRES )
			createDomainsMultiResolution();

		fprintf( stderr, "%s: %d: %d domains created ...\n", __FILE__, __LINE__, nTotalDomain );
	}

	void createDomainsMultiSize()
	{
		float lowF[nDim];
		float highF[nDim];
		int dim[nDim], dsfDim[nDim];
		int resampledFieldSize[4];
		T* resampledData = NULL;

		fprintf( stderr, "Creating domains of different sizes from single resolution ... \n");

		// Down sample data to one level
		fprintf( stderr, "Down sampling field ... \n");
		downsampledDataField = new ITL_field_regular<T>();
		downsampledBinField = new ITL_field_regular<int>();

		//inputDataField->getSize( dim );
		(*inputDataField)->getSize( dim );
		FEL_util<T>::resampleField( (*inputDataField)->getDataFull(), dataType, dim,
									&resampledData, resampledFieldSize,
									0.5f );

		// Create a down sampled data field
		fprintf( stderr, "Creating a field class from down sampled data ... \n");
		assert( resampledData != NULL );
		for( int i=0; i<nDim; i++ )
		{
			lowF[i] = 0;
			highF[i] = resampledFieldSize[i]-1;
		}
		downsampledDataField->initialize( resampledData, nDim, lowF, highF );

		// Create a down sampled bin field
		fprintf( stderr, "Computing histogram bins from down sampled data ... \n");
		downsampledBinField->initialize( nDim, lowF, highF );
		downsampledBinField->getSize( dsfDim );
		if( histMappingType == 2 )
			FEL_util<T>::computeHistogramBinField2D( dataType,
													 &downsampledDataField,
													 &downsampledBinField,
													 histogram, nBin,
													 dsfDim );
		else
			FEL_util<T>::computeHistogramBinField( dataType,
												   &downsampledDataField,
												   &downsampledBinField,
												   histogram, nBin,
												   dsfDim );

		// Create regular size domains
		// of three different sizes
		// BxB, 2Bx2B, 4Bx4B
		// Take special care for 2D fields
		fprintf( stderr, "Creating domains of multiple sizes ... \n");
		int domainSize[4];
		for( int i=0; i<3; i++ )
		{
			// Set domain size for this iteration
			domainSize[0] = smallestDomainSize[0] * pow( 2, i);
			domainSize[1] = smallestDomainSize[1] * pow( 2, i);
			if( dsfDim[2] == 2 )
				domainSize[2] = smallestDomainSize[2];
			else
				domainSize[2] = smallestDomainSize[2] * pow( 2, i);

			// Create domains
			fprintf( stderr, "%s: %d: Creating domains of size <%d, %d, %d> ...\n", __FILE__, __LINE__,
							 domainSize[0], domainSize[1], domainSize[2] );
			createRegularDomains( domainSize, i );

		}// end for : i

		// Clear up
		delete downsampledDataField;
		delete downsampledBinField;
	}

	void createDomainsMultiResolution()
	{
		float lowF[nDim];
		float highF[nDim];
		int dim[nDim], dsfDim[nDim];
		int resampledFieldSize[3];
		T* resampledData = NULL;

		fprintf( stderr, "Creating domains of same sizes from different resolutions ... \n");
		clock_t starttime = ITL_util<float>::startTimer();

		// Create regular size domains
		// from fields of three different resolutions
		// (N/2)^3, (N/4)^3 and (N/8)^3
		downsampledDataFieldArray = new ITL_field_regular<T>[3];
		downsampledBinFieldArray = new ITL_field_regular<int>[3];

		for( int i=0; i<3; i++ )
		{
			// Down-sample the field
			fprintf( stderr, "%d: Down sampling vector field ... \n", i );
			//if( resampledData != NULL ) delete resampledData;

			if( i == 0 )
			{
				//inputDataField->getSize( dim );
				(*inputDataField)->getSize( dim );
				FEL_util<T>::resampleField( (*inputDataField)->getDataFull(), dataType, dim,
											&resampledData, resampledFieldSize,
											0.5f );
			}
			else
			{
				downsampledDataFieldArray[i-1].getSize( dim );
				FEL_util<T>::resampleField( downsampledDataFieldArray[i-1].getDataFull(), dataType, dim,
											&resampledData, resampledFieldSize,
											0.5f );
			}

			fprintf( stderr, "%d: Down sampling done ... \n", i );
			// Copy down-sampled data
			// and create a field
			int nPointField = 1;
			for( int j=0; j<nDim; j++ )
			{
				nPointField *= resampledFieldSize[j];
				lowF[j] = 0;
				highF[j] = resampledFieldSize[j]-1;
			}

			// Create a down sampled data field
			fprintf( stderr, "%d: Creating a field class from down sampled data ... \n", i );
			downsampledDataFieldArray[i].initialize( resampledData, nDim, lowF, highF );
			delete [] resampledData;

			// Create a down sampled bin field
			fprintf( stderr, "%d: Computing histogram bins from down sampled data ... \n", i );
			downsampledBinFieldArray[i].initialize( nDim, lowF, highF );
			downsampledBinFieldArray[i].getSize( dsfDim );
			ITL_field_regular<T>* tmp = downsampledDataFieldArray+i;
			ITL_field_regular<int>* tmp2 = downsampledBinFieldArray+i;
			if( histMappingType == 2 )
				FEL_util<T>::computeHistogramBinField2D( dataType,
														 &tmp,
														 &tmp2,
														 histogram, nBin,
														 dsfDim );
			else
				FEL_util<T>::computeHistogramBinField( dataType,
													   &tmp,
													   &tmp2,
													   histogram, nBin,
													   dsfDim );

			// Create domains
			fprintf( stderr, "Creating domains of size <%d, %d, %d> ...\n", smallestDomainSize[0], smallestDomainSize[1], smallestDomainSize[2] );
			createRegularDomains2( smallestDomainSize, i );

		}// end for : i

		// Create additional domains through adaptive partitioning
		// *********To be incorporated************

		domainCreationTime = ITL_util<float>::endTimer( starttime );
	}

	void createRegularDomains( int* domainSize, int sizeIndex )
	{
		int nDomain[nDim];
		float lowSub[nDim];
		float highSub[nDim];

		// Compute the number of domains
		int dsfDim[nDim];
		downsampledDataField->getSize( dsfDim );
		nDomain[0] = (int)ceil( dsfDim[0]/(float)domainSize[0] );
		nDomain[1] = (int)ceil( dsfDim[1]/(float)domainSize[1] );
		nDomain[2] = (int)ceil( dsfDim[2]/(float)domainSize[2] );

		#ifdef DEBUG_MODE
		printf( "Dimension of the downsampled field: %d, %d, %d\n", dsfDim[0], dsfDim[1], dsfDim[2] );
		printf( "Number of domains to be created: %d, %d, %d\n", nDomain[0], nDomain[1], nDomain[2] );
		printf( "Size of a domain: %d, %d, %d\n", domainSize[0], domainSize[1], domainSize[2] );
		#endif

		// Nested loop for partitioning
		int dsfLowInt[nDim];
		int dsfHighInt[nDim];
		downsampledDataField->getBounds( dsfLowInt, dsfHighInt );

		int domainIndex = nTotalDomain;

		for( int z=0; z<nDomain[2]; z++ )
		{
			// Determine [zmin, zmax] for the next block
			if( dsfDim[2] == 2 )
			{
				lowSub[2] = 0; highSub[2] = 1;
			}
			else
			{
				lowSub[2] = z * domainSize[2];
				if( isSharingNeighbor == true )
					highSub[2] = std::min( (float)dsfHighInt[2], lowSub[2] + domainSize[2] );
				else
					highSub[2] = std::min( (float)dsfHighInt[2], lowSub[2] + domainSize[2] - 1 );
			}

			for(int y=0; y<nDomain[1]; y++ )
			{
				// Determine [ymin, ymax] for the next block
				lowSub[1] = y * domainSize[1];
				if( isSharingNeighbor == true )
					highSub[1] = std::min( (float)dsfHighInt[1], lowSub[1] + domainSize[1] );
				else
					highSub[1] = std::min( (float)dsfHighInt[1], lowSub[1] + domainSize[1] - 1 );

				for( int x=0; x<nDomain[0]; x++)
				{
					// Determine [xmin, xmax] for the next block
					lowSub[0] = x * domainSize[0];
					if( isSharingNeighbor == true )
						highSub[0] = std::min( (float)dsfHighInt[0], lowSub[0] + domainSize[0] );
					else
						highSub[0] = std::min( (float)dsfHighInt[0], lowSub[0] + domainSize[0] - 1 );

					// clamp domain limits within field
					//lowSub[0] = ITL_util<float>::clamp( 0, dsfDim[0]-1 );
					//highSub[0] = ITL_util<float>::clamp( 0, dsfDim[0]-1 );
					//lowSub[1] = ITL_util<float>::clamp( 0, dsfDim[1]-1 );
					//highSub[1] = ITL_util<float>::clamp( 0, dsfDim[1]-1 );
					//lowSub[2] = ITL_util<float>::clamp( 0, dsfDim[2]-1 );
					//highSub[2] = ITL_util<float>::clamp( 0, dsfDim[2]-1 );


					// Set properties of next domain block
					//domainList[domainIndex].initialize( nDim, lowSub, highSub );
					FEL_domain<T> newDomain( nDim, lowSub, highSub, sizeIndex, nBin );

					// Get domain block size
					int nPointDomain = newDomain.getSize();

					// Get domain block data
					int* domainBinData = new int[nPointDomain];
					downsampledBinField->getDataBetween( lowSub, highSub, domainBinData );
					#ifdef DEBUG_MODE
					int m = ITL_util<int>::Min( domainBinData, nPointDomain );
					int M = ITL_util<int>::Max( domainBinData, nPointDomain );
					assert( m >= 0 );
					assert( M < nBin );
					#endif

					// Convert to distribution
					double* domainFreqList = new double[nBin];
					FEL_util<double>::computeHistogramFromBinField( domainBinData, nPointDomain, nBin, domainFreqList );

					// Store computed domain distribution
					newDomain.setDistribution( domainFreqList, nBin );
					#ifdef DEBUG_MODE
					double m = ITL_util<double>::Min( domainFreqList, nBin );
					double M = ITL_util<double>::Max( domainFreqList, nBin );
					fprintf( stderr, "Domain distribution min/max: %g, %g\n", m, M );
					#endif

					// Store domain in list
					domainList.push_back( newDomain );

					// Store domain extent for matlab plotting
					//domainLimitList[domainIndex*7] = lowSub[0];
					//domainLimitList[domainIndex*7+1] = highSub[0];
					//domainLimitList[domainIndex*7+2] = lowSub[1];
					//domainLimitList[domainIndex*7+3] = highSub[1];
					//domainLimitList[domainIndex*7+4] = lowSub[2];
					//domainLimitList[domainIndex*7+5] = highSub[2];
					//domainLimitList[domainIndex*7+6] = 1;

					// Increment the index to the next domain
					domainIndex ++;

					// Free temporary resources
					delete [] domainBinData;
					delete [] domainFreqList;
				}
			}
		}

		// Update total number of domains
		nTotalDomain += ITL_util<int>::prod( nDomain, nDim );
	}

	void
	createRegularDomains2( int* domainSize, int lodIndex )
	{
		int nDomain[nDim];
		float lowSub[nDim];
		float highSub[nDim];

		// Compute the number of domains
		int dsfDim[nDim];
		downsampledBinFieldArray[lodIndex].getSize( dsfDim );
		nDomain[0] = (int)ceil( dsfDim[0]/(float)domainSize[0] );
		nDomain[1] = (int)ceil( dsfDim[1]/(float)domainSize[1] );
		nDomain[2] = (int)ceil( dsfDim[2]/(float)domainSize[2] );

		//#ifdef DEBUG_MODE
		printf( "Dimension of the field: %d, %d, %d\n", dsfDim[0], dsfDim[1], dsfDim[2] );
		printf( "Number of domains to be created: %d, %d, %d\n", nDomain[0], nDomain[1], nDomain[2] );
		printf( "Size of a domain: %d, %d, %d\n", domainSize[0], domainSize[1], domainSize[2] );
		//#endif

		// Nested loop for partitioning
		int dsfLowInt[nDim];
		int dsfHighInt[nDim];
		downsampledBinFieldArray[lodIndex].getBounds( dsfLowInt, dsfHighInt );
		//#ifdef DEBUG_MODE
		printf( "Extent of downsampled field: <%d,%d,%d>-<%d,%d,%d>\n", dsfLowInt[0], dsfLowInt[1], dsfLowInt[2],
																		dsfHighInt[0], dsfHighInt[1], dsfHighInt[2] );
		//#endif

		int domainIndex = nTotalDomain;

		for( int z=0; z<nDomain[2]; z++ )
		{
			// Determine [zmin, zmax] for the next block
			if( nDomain[2] == 2 )
			{
				lowSub[2] = 0; highSub[2] = 1;
			}
			else
			{
				lowSub[2] = z * domainSize[2];
				if( isSharingNeighbor == true )
					highSub[2] = std::min( (float)dsfHighInt[2], lowSub[2] + domainSize[2] );
				else
					highSub[2] = std::min( (float)dsfHighInt[2], lowSub[2] + domainSize[2] - 1 );
			}

			for(int y=0; y<nDomain[1]; y++ )
			{
				// Determine [ymin, ymax] for the next block
				lowSub[1] = y * domainSize[1];
				if( isSharingNeighbor == true )
					highSub[1] = std::min( (float)dsfHighInt[1], lowSub[1] + domainSize[1] );
				else
					highSub[1] = std::min( (float)dsfHighInt[1], lowSub[1] + domainSize[1] - 1 );

				for( int x=0; x<nDomain[0]; x++)
				{
					// Determine [xmin, xmax] for the next block
					lowSub[0] = x * domainSize[0];
					if( isSharingNeighbor == true )
						highSub[0] = std::min( (float)dsfHighInt[0], lowSub[0] + domainSize[0] );
					else
						highSub[0] = std::min( (float)dsfHighInt[0], lowSub[0] + domainSize[0] - 1 );

					// clamp domain limits within field
					//lowSub[0] = ITL_util<float>::clamp( 0, dsfDim[0]-1 );
					//highSub[0] = ITL_util<float>::clamp( 0, dsfDim[0]-1 );
					//lowSub[1] = ITL_util<float>::clamp( 0, dsfDim[1]-1 );
					//highSub[1] = ITL_util<float>::clamp( 0, dsfDim[1]-1 );
					//lowSub[2] = ITL_util<float>::clamp( 0, dsfDim[2]-1 );
					//highSub[2] = ITL_util<float>::clamp( 0, dsfDim[2]-1 );

					// Set properties of next domain block
					//domainList[domainIndex].initialize( nDim, lowSub, highSub );
					FEL_domain<T> newDomain( nDim, lowSub, highSub, lodIndex, nBin );

					// Get domain block size
					int nPointDomain = newDomain.getSize();

					// Get domain bin data
					int* domainBinData = new int[nPointDomain];
					downsampledBinFieldArray[lodIndex].getDataBetween( lowSub, highSub, domainBinData );

					// Convert to distribution
					double* domainFreqList = new double[nBin];
					FEL_util<double>::computeHistogramFromBinField( domainBinData, nPointDomain, nBin, domainFreqList );

					// Store computed domain distribution
					newDomain.setDistribution( domainFreqList, nBin );

					// Store domain in list
					domainList.push_back( newDomain );

					// Store domain extent for matlab plotting
					//domainLimitList[domainIndex*7] = lowSub[0];
					//domainLimitList[domainIndex*7+1] = highSub[0];
					//domainLimitList[domainIndex*7+2] = lowSub[1];
					//domainLimitList[domainIndex*7+3] = highSub[1];
					//domainLimitList[domainIndex*7+4] = lowSub[2];
					//domainLimitList[domainIndex*7+5] = highSub[2];
					//domainLimitList[domainIndex*7+6] = 1;

					// Increment the index to the next domain
					domainIndex ++;

					// Free temporary resources
					delete [] domainBinData;
				}
			}
		}

		// Update total number of domains
		nTotalDomain += ITL_util<int>::prod( nDomain, nDim );
	}

	void
	computeDomainEntropies()
	{
		int domainBounds[6];
		int nPointDomain = 1;
		double domainFreqList[nBin];
		float domainEntropy;
		int nb;

		int domainID = 0;
		for( typename list< FEL_domain<T> >::iterator domainIter = domainList.begin();
			 domainIter != domainList.end();
			 domainIter++ )
		{
			cout << domainID << endl;
			domainIter->getBounds( domainBounds );

			cout << "1" << endl;
			domainIter->getDistribution( domainFreqList, &nb );

			// Compute and store domain entropy
			cout << "2" << endl;
			domainEntropy = (float)ITL_entropycore::computeEntropy_HistogramBased2( domainFreqList, nBin, false );
			domainIter->setEntropy( domainEntropy );

			if( domainID == 0 )
			{
				minDomainEntropy = maxDomainEntropy = domainEntropy;
			}
			else
			{
				if( domainEntropy < minDomainEntropy ) minDomainEntropy = domainEntropy;
				if( domainEntropy > maxDomainEntropy ) maxDomainEntropy = domainEntropy;
			}

			cout << "3" << endl;
			// Store domain entropies for display
			if( matlab_plotting_on )
			{
				plotter->domainEntropyList[domainID] = domainEntropy;
			}


			domainID ++;

			// Free temporary resources
			//delete [] domainBinData;
		}
	}

	void
	sortDomains_Entropybased()
	{
		// Sort domains based on entropy
		domainList.sort( FEL_domain<T>::compare_entropy );
	}

	void
	pruneDomains_Entropybased()
	{
		int i = 0;
		list< FEL_domain<T> > reducedDomainList;

		float domainEntropyRange = maxDomainEntropy - minDomainEntropy;
		float entropyThreshold = 0.01f * domainEntropyRange;

		float lastEntropy, curEntropy;
		for( typename list< FEL_domain<T> >::iterator domainIter = domainList.begin();
			 domainIter != domainList.end();
			 domainIter++ )
		{
			// Get entropy of current domain
			curEntropy = domainIter->getEntropy();

			if( i==0 )
			{
				lastEntropy = curEntropy;
				reducedDomainList.push_back( *domainIter );
			}
			else
			{
				if( ( curEntropy - lastEntropy ) > entropyThreshold )
				{
					lastEntropy = curEntropy;
					reducedDomainList.push_back( *domainIter );
				}
			}

			i++;

		}// end for

		// Update list of domains
		domainList.clear();
		domainList = reducedDomainList;

	}

	void
	pruneDomains_Entropybased2()
	{
		int i = 0;
		list< FEL_domain<T> > curDomainList = domainList;
		list< FEL_domain<T> > reducedDomainList;

		// Determine the pruning thresholds
		float domainEntropyRange = maxDomainEntropy - minDomainEntropy;
		float entropyThreshold = 0.01f * domainEntropyRange;
		double distanceThreshold = 0.01f;

		float lastEntropy, curEntropy;
		double lastDistribution[nBin], curDistribution[nBin];

		// Iterate over the domain set until each distribution is different from its neighbors
		int nStartDomain = curDomainList.size();
		int iter = 0;
		int iDomain = 0;
		int nReduced = 0;
		int nb;

		do
		{
			// Clear the reduced domain list
			reducedDomainList.clear();

			// Reset indices
			iDomain = 0;
			nReduced = 0;

			// Next Iteration
			for( typename list< FEL_domain<T> >::iterator domainIter = curDomainList.begin();
				 domainIter != curDomainList.end();
				 domainIter++ )
			{
				// Get distribution and entropy of current domain
				domainIter->getDistribution( curDistribution, &nb );
				curEntropy = domainIter->getEntropy();

				if( iDomain == 0 )
				{
					lastEntropy = curEntropy;
					memcpy( lastDistribution, curDistribution, sizeof(double)*nBin );
					reducedDomainList.push_back( *domainIter );
				}
				else
				{
					if( ( curEntropy - lastEntropy ) > entropyThreshold )
					{
						// Keep domain if distributions have too different entropy level
						lastEntropy = curEntropy;
						memcpy( lastDistribution, curDistribution, sizeof(double)*nBin );
						reducedDomainList.push_back( *domainIter );
					}
					else if( FEL_util<double>::computeHistogramMatchError( distanceMeasureType,
																		   lastDistribution, curDistribution, nBin ) > distanceThreshold )
					{
						// Keep domain if entropies are similar, but distributions are different
						lastEntropy = curEntropy;
						memcpy( lastDistribution, curDistribution, sizeof(double)*nBin );
						reducedDomainList.push_back( *domainIter );
					}
					else
					{
						nReduced ++;
					}
				}// end if-else : first domain or not

				// Go to next domain
				iDomain++;

			}// end for

			// The reduced set of domain becomes current for next iteration
			curDomainList.clear();
			curDomainList = reducedDomainList;

			fprintf( stderr, "iter-%d: Number of discarded/remaining domains: %d/%d\n", iter, nReduced, reducedDomainList.size() );
			iter++;

		}while( nReduced > 0 );// end do-while

		// Update list of domains
		domainList.clear();
		domainList = reducedDomainList;

	}


	void organizeDomains()
	{
		clock_t starttime = ITL_util<float>::startTimer();

		// Compute domain entropies
		fprintf( stderr, "%s: %d: Computing domain entropies ...\n", __FILE__, __LINE__ );
		computeDomainEntropies();

		// Sort domains based on entropy
		fprintf( stderr, "%s: %d: Sorting domains based on entropy ...\n", __FILE__, __LINE__ );
		sortDomains_Entropybased();
		fprintf( stderr, "%s: %d: Domain entropy range: %g %g ...\n", __FILE__, __LINE__, minDomainEntropy, maxDomainEntropy );

		// Prune redundant domains
		if( isDomainPruningOn == true )
		{
			fprintf( stderr, "%s: %d: Reducing domain set based on entropy ...\n", __FILE__, __LINE__ );
			fprintf( stderr, "Number of domains before pruning: %d\n", domainList.size() );
			//pruneDomains_Entropybased();
			pruneDomains_Entropybased2();
			fprintf( stderr, "Number of domains after pruning: %d\n", domainList.size() );

			// Update the total number of domains
			nTotalDomain = domainList.size();
		}

		// Organize sorted domains in a segment tree
		fprintf( stderr, "%s: %d: Storing domains in a hierarchy ...\n", __FILE__, __LINE__ );
		domainTree.initTree( domainList );
		domainTree.createTree( domainTree.getRoot(), 0, nTotalDomain-1, 0 );

		domainOrgTime = ITL_util<float>::endTimer( starttime );
		fprintf( stderr, "Domain Organization time: %g\n", domainOrgTime );
	}

	FEL_segmenttree<T>*
	getDomainTree()
	{
		return (&domainTree);
	}

	void
	encode()
	{
		clock_t starttime, endtime;
		float lowF[nDim], highF[nDim];
		int dim[nDim];

		if( encodingType == REGULAR )
		{
			fprintf( stderr, "Initializing encoding with regular partitioning ... \n" );
			starttime = ITL_util<float>::startTimer();
			encode_Regular();
			encodeTime = ITL_util<float>::endTimer( starttime );
		}
		else if( encodingType == QUADTREE )
		{
			fprintf( stderr, "Initializing encoding with quadtree partitioning ... \n" );

			inputBinField->getBounds( lowF, highF );
			inputBinField->getSize( dim );
			starttime = ITL_util<float>::startTimer();
			encode_Quadtree2( lowF, highF, dim );
			encodeTime = ITL_util<float>::endTimer( starttime );

		}
		else if( encodingType == OCTREE )
		{
			fprintf( stderr, "Initializing encoding with octree partitioning ... \n" );

			inputBinField->getBounds( lowF, highF );
			inputBinField->getSize( dim );
			starttime = ITL_util<float>::startTimer();
			encode_Octree2( lowF, highF, dim );
			encodeTime = ITL_util<float>::endTimer( starttime );
		}

		// Analysis
		double minError = 100;
		double sumError = 0;
		double avgError = 0;
		double maxError = -100;
		for( int iR=0; iR<nLeafRange; iR++ )
		{
			if( plotter->errorList[iR] < minError )	minError = plotter->errorList[iR];
			if( plotter->errorList[iR] > maxError )	maxError = plotter->errorList[iR];
			sumError += plotter->errorList[iR];
		}
		avgError = ( sumError / (double)nLeafRange );

		double minMuError = 100;
		double sumMuError = 0;
		double avgMuError = 0;
		double maxMuError = -100;
		for( int iR=0; iR<nLeafRange; iR++ )
		{
			if( plotter->muErrorList[iR] < minMuError )	minMuError = plotter->muErrorList[iR];
			if( plotter->muErrorList[iR] > maxMuError )	maxMuError = plotter->muErrorList[iR];
			sumMuError += plotter->muErrorList[iR];
		}
		avgMuError = ( sumMuError / (double)nLeafRange );

		double minSdError = 100;
		double sumSdError = 0;
		double avgSdError = 0;
		double maxSdError = -100;
		for( int iR=0; iR<nLeafRange; iR++ )
		{
			if( plotter->sdErrorList[iR] < minSdError )	minSdError = plotter->sdErrorList[iR];
			if( plotter->sdErrorList[iR] > maxSdError )	maxSdError = plotter->sdErrorList[iR];
			sumSdError += plotter->sdErrorList[iR];
		}
		avgSdError = ( sumSdError / (double)nLeafRange );

		fprintf( stderr, "Total encoding time: %g\n", encodeTime );
		fprintf( stderr, "Min/avg/max match error: %g %g %g\n", minError, avgError, maxError );
		fprintf( stderr, "Min/avg/max Mean difference error: %g %g %g\n", minMuError, avgMuError, maxMuError );
		fprintf( stderr, "Min/avg/max SD difference error: %g %g %g\n", minSdError, avgSdError, maxSdError );

	}// End function

	void
	encode_Regular()
	{
		int nRange[nDim];
		float lowSub[nDim];
		float highSub[nDim];
		int lowInt[nDim], highInt[nDim];
		int dim[nDim];
		inputBinField->getBounds( lowInt, highInt );
		inputBinField->getSize( dim );

		// Compute the number of ranges
		nRange[0] = (int)ceil( dim[0]/(float)smallestRangeSize[0] );
		nRange[1] = (int)ceil( dim[1]/(float)smallestRangeSize[1] );
		nRange[2] = (int)ceil( dim[2]/(float)smallestRangeSize[2] );

		//#ifdef DEBUG_MODE
		printf( "Dimension of the field: %d, %d, %d\n", dim[0], dim[1], dim[2] );
		printf( "Number of ranges to be created: %d, %d, %d\n", nRange[0], nRange[1], nRange[2] );
		printf( "Size of a range: %d, %d, %d\n", smallestRangeSize[0], smallestRangeSize[1], smallestRangeSize[2] );
		//#endif

		// Nested loop for partitioning
		// data into ranges for t-th time step
		int rangeIndex = 0;
		int nPointRange = 0;
		double rangeFreqList[nBin];
		for( int z=0; z<nRange[2]; z++ )
		{
			// Determine [zmin, zmax] for the next range
			if( nRange[2] == 2 )
			{
				lowSub[2] = 0; highSub[2] = 1;
			}
			else
			{
				lowSub[2] = z * smallestRangeSize[2];
				if( isSharingNeighbor == true )
					highSub[2] = std::min( (float)highInt[2], lowSub[2] + smallestRangeSize[2] );
				else
					highSub[2] = std::min( (float)highInt[2], lowSub[2] + smallestRangeSize[2]-1 );
			}

			for(int y=0; y<nRange[1]; y++ )
			{
				// Determine [ymin, ymax] for the next block
				lowSub[1] = y * smallestRangeSize[1];
				if( isSharingNeighbor == true )
					highSub[1] = std::min( (float)highInt[1], lowSub[1] + smallestRangeSize[1] );
				else
					highSub[1] = std::min( (float)highInt[1], lowSub[1] + smallestRangeSize[1]-1 );

				for( int x=0; x<nRange[0]; x++)
				{
					// Determine [xmin, xmax] for the range
					lowSub[0] = x * smallestRangeSize[0];
					if( isSharingNeighbor == true )
						highSub[0] = std::min( (float)highInt[0], lowSub[0] + smallestRangeSize[0] );
					else
						highSub[0] = std::min( (float)highInt[0], lowSub[0] + smallestRangeSize[0]-1 );

					// Get range block data and convert it to distribution
					nPointRange = 1;
					for( int i=0; i<nDim; i++ )
						nPointRange *= (int)( highSub[i] - lowSub[i] + 1 );
					int* rangeBinData = new int[nPointRange];
					inputBinField->getDataBetween( lowSub, highSub, rangeBinData );
					FEL_util<double>::computeHistogramFromBinField( rangeBinData, nPointRange, nBin, rangeFreqList );
					delete [] rangeBinData;

					// Create a range block and
					// set properties of the range block including distribution
					FEL_range<T> nextRange( nDim, lowSub, highSub, nBin );
					nextRange.setDistribution( rangeFreqList, nBin );

					// Add new range to the list
					rangeBlockList.push_back( nextRange );

					// Store domain extent for matlab plotting
					if( matlab_plotting_on )
					{
						plotter->rangeLimitList[rangeIndex*7] = lowSub[0];
						plotter->rangeLimitList[rangeIndex*7+1] = highSub[0];
						plotter->rangeLimitList[rangeIndex*7+2] = lowSub[1];
						plotter->rangeLimitList[rangeIndex*7+3] = highSub[1];
						plotter->rangeLimitList[rangeIndex*7+4] = lowSub[2];
						plotter->rangeLimitList[rangeIndex*7+5] = highSub[2];
						plotter->rangeLimitList[rangeIndex*7+6] = 1;
					}

					// Increment the index to the next domain
					rangeIndex ++;
				}
			}
		}

		// Update total number of ranges
		nTotalRange = ITL_util<int>::prod( nRange, nDim );
		fprintf( stderr, "%s: %d: %d ranges created ...\n", __FILE__, __LINE__, nTotalRange );
		fprintf( stderr, "%s: %d: starting encoder ...\n", __FILE__, __LINE__ );
		//printRanges();

		// Iterate though each range
		// and find the nearest domain
		int bestMatch = -1;
		double matchError = -1;
		double muError = -1;
		double sdError = -1;

		int rangeID = 0;
		int optRot = 0;
		bool optRef = false;
		float rangeLow[3];
		float rangeHigh[3];

		// Iterate through the ranges
		// For each range, scan through ALL
		// or similar entropy domains
		for( typename list< FEL_range<T> >::iterator rangeBlockIter = rangeBlockList.begin();
			 rangeBlockIter != rangeBlockList.end();
			 rangeBlockIter++ )
		{
			fprintf( stderr, "Finding best match for the %d-th next range ... \n", rangeID );

			// Get Range distribution
			rangeBlockIter->getDistribution( rangeFreqList, &nBin );
			#ifdef DEBUG_MODE
			double m = ITL_util<double>::Min( rangeFreqList, nBin );
			double M = ITL_util<double>::Max( rangeFreqList, nBin );
			fprintf( stderr, "%d-th Range distribution min/max: %g, %g\n", rangeID, m, M );
			#endif

			//bestMatch = matchRangeToDomains( rangeBlockIter, &matchError, &optimalTrMatrix );
			optRot = 0;
			optRef = false;
			bestMatch = matchRangeToDomains2( rangeFreqList, &matchError, &optRot, &optRef, &muError, &sdError );

			#ifdef DEBUG_MODE
			fprintf( stderr, "%d-th range maps to %d-th domain with %g <%d, %d> error ...\n",
					 rangeID, bestMatch, matchError, optRot, (int)optRef );
			#endif

			// Add range information to domain
			fprintf( stderr, "Adding range information to the matched domain ... \n" );
			typename list<FEL_domain<T> >::iterator tmpIter = domainList.begin();
			for( int iD=0; iD<bestMatch; iD++ )
				tmpIter++;
			tmpIter->addMatchingRange( rangeID, matchError, optRot, optRef );

			// Add domain information to range
			fprintf( stderr, "Adding domain information to the matched range ... \n" );
			rangeBlockIter->setMatchedDomainId( bestMatch );
			rangeBlockIter->setMatchError( matchError  );
			rangeBlockIter->setMatchTransformation( optRot, optRef );

			// Store Matching information
			// for Matlab plotting
			if( matlab_plotting_on )
			{
				// Record matching error
				plotter->errorList[rangeID] = matchError;
				plotter->muErrorList[rangeID] = muError;
				plotter->sdErrorList[rangeID] = sdError;

				// Record entropy of each range and it matching domain
				plotter->rangeEntropyList[rangeID] = (float)ITL_entropycore::computeEntropy_HistogramBased2( rangeFreqList, nBin, false );
				plotter->matchedDomainIdList[rangeID] = bestMatch;

				plotter->domainMatchFreqList[bestMatch] ++;
			}

			// Increment range ID
			rangeID ++;

		}// end for : scan through ranges

		if( matlab_plotting_on )
		{
			plotter->nBin = nBin;
			plotter->nDomain = domainList.size();
			plotter->nRange = rangeBlockList.size();
		}


		for( typename list< FEL_domain<T> >::iterator domIter = domainList.begin();
			 domIter != domainList.end();
			 domIter++ )
		{
			domIter->sortCodeBookEntries();
		}

	}

	void
	encode_Quadtree2( float *parentLow, float *parentHigh, int *rangeSize  )
	{
		// Initialize range block tree
		// Send pointers to related variables
		float lowData[nDim];
		float highData[nDim];

		//inputDataField->getBounds( lowData, highData );
		(*inputDataField)->getBounds( lowData, highData );

		float lowBin[nDim];
		float highBin[nDim];
		inputBinField->getBounds( lowBin, highBin );
		ITL_field_regular<int>* tmp2 = new ITL_field_regular<int>( inputBinField->getDataFull(),
																   nDim, lowBin, highBin );

		rangeBlockTree.initialize( &domainList,
							 	   &histogram, nBin,
								   &tmp2,
								   smallestRangeSize,
								   distanceMeasureType,
								   encodingType,
								   matchErrorThreshold,
								   entropyThreshold,
								   &trManager );

		rangeBlockTree.setMatlabPlotter( &plotter, matlab_plotting_on );

		// Create root node (represents entire field)
		fprintf( stderr, "Creating range blobk quadtree root ... \n");
		rangeBlockTree.initTree( lowData, highData );

		// Recursive create range blocks
		// Find matching domains and store ranges in tree
		fprintf( stderr, "Recursively creating range block quadtree ... \n");
		rangeBlockTree.createTree2( rangeBlockTree.getRoot(), lowData, highData, 0 );

		// Obtain total number of ranges
		nTotalRange = rangeBlockTree.getNumRangeBlock();
		nLeafRange = rangeBlockTree.getNumLeafRangeBlock();

		if( matlab_plotting_on )
		{
			plotter->nBin = nBin;
			plotter->nDomain = domainList.size();
			plotter->nRange = nLeafRange;
		}

		fprintf( stderr, "%s: %d: %d rangeblocks created ...\n", __FILE__, __LINE__, nTotalRange );
		fprintf( stderr, "%s: %d: %d leaf nodes created ...\n", __FILE__, __LINE__, nLeafRange );
		fprintf( stderr, "%s: %d: No good match for %d ranges ...\n", __FILE__, __LINE__, rangeBlockTree.noGoodMatchFound );

		//delete tmp;
		delete tmp2;
	}

	void
	encode_Octree2( float *parentLow, float *parentHigh, int *rangeSize  )
	{
		float lowData[nDim];
		float highData[nDim];
		//inputDataField->getBounds( lowData, highData );
		(*inputDataField)->getBounds( lowData, highData );

		float lowBin[nDim];
		float highBin[nDim];
		inputBinField->getBounds( lowBin, highBin );
		ITL_field_regular<int>* tmp2 = new ITL_field_regular<int>( inputBinField->getDataFull(),
																   nDim, lowBin, highBin );

		// Initialize range block tree
		// Send pointers to related variables
		rangeBlockTree.initialize( &domainList,
								   &histogram, nBin,
								   &tmp2,
								   //&inputDataField, &inputBinField,
								   smallestRangeSize,
								   distanceMeasureType,
								   encodingType,
								   matchErrorThreshold,
								   entropyThreshold,
								   &trManager );

		rangeBlockTree.setMatlabPlotter( &plotter, matlab_plotting_on );

		// Create root node (represents entire field)
		fprintf( stderr, "Creating range blobk quadtree root ... \n");
		rangeBlockTree.initTree( lowData, highData );

		// Recursive create range blocks
		// Find matching domains and store ranges in tree
		fprintf( stderr, "Recursively creating range block octree ... \n");
		rangeBlockTree.createTree2( rangeBlockTree.getRoot(), lowData, highData, 0 );

		// Obtain total number of ranges
		nTotalRange = rangeBlockTree.getNumRangeBlock();
		nLeafRange = rangeBlockTree.getNumLeafRangeBlock();

		if( matlab_plotting_on )
		{
			plotter->nBin = nBin;
			plotter->nDomain = domainList.size();
			plotter->nRange = nLeafRange;
			plotter->nUnmatchedRange = rangeBlockTree.noGoodMatchFound;

		}

		int nSavedDomain = 0;
		for( typename list< FEL_domain<T> >::iterator domainIter = domainList.begin();
			 domainIter != domainList.end();
			 domainIter++ )
		{
			if( domainIter->getNumMatchingRange() > 0 )
				nSavedDomain ++;
		}

		fprintf( stderr, "%s: %d: %d domainblocks utilized ...\n", __FILE__, __LINE__, nSavedDomain );
		fprintf( stderr, "%s: %d: %d rangeblocks created ...\n", __FILE__, __LINE__, nTotalRange );
		fprintf( stderr, "%s: %d: %d leaf nodes created ...\n", __FILE__, __LINE__, nLeafRange );
		fprintf( stderr, "%s: %d: No good match for %d ranges ...\n", __FILE__, __LINE__, rangeBlockTree.noGoodMatchFound );
		double spaceSavingFromData = 1 -  ( (double)( (rangeBlockTree.noGoodMatchFound + nSavedDomain)*nBin*4 + nSavedDomain*5 ) /
											(double)( (*inputDataField)->getSize()*4 )
										  );

		double spaceSavingFromHist = 1 -  ( (double)( (rangeBlockTree.noGoodMatchFound + nSavedDomain)*nBin*4 + nSavedDomain*5 ) /
											(double)( nTotalRange*nBin*4 )
										  );
		fprintf( stderr, "%s: %d: Space savings: %g %g\n", __FILE__, __LINE__, spaceSavingFromData,spaceSavingFromHist );
	}

	void
	encode_Regular_Distribution( double** rangeHistList, int nRange )
	{
		int nb;
		double* rangeBinData = new double[nBin];
		int bestMatch;
		int optRot;
		bool optRef;
		double matchError, muError, sdError;
		float lowSub[3] = {0,0,0};
		float highSub[3] = {3, 3, 3};

		//#ifdef DEBUG_MODE
		fprintf( stderr, "Number of domain and range histograms: %d, %d\n", nTotalDomain, nRange );
		//#endif

		// Put all ranges in a list
		for( int iR=0; iR<nRange; iR++ )
		{
			// Copy range distribution
			memcpy( rangeBinData, rangeHistList[iR], sizeof(double)*nBin );

			// Create a range block and
			// set properties of the range block including distribution
			FEL_range<T> nextRange( 3, lowSub, highSub, nBin );
			nextRange.setDistribution( rangeBinData, nBin );

			// Add new range to the list
			rangeBlockList.push_back( nextRange );
		}

		// Loop for encoding
		typename list< FEL_range<T> >::iterator rangeBlockIter = rangeBlockList.begin();
		for( int iR=0; iR<nRange; iR++ )
		{
			#ifdef DEBUG_MODE
			fprintf( stderr, "Finding best match for the %d-th range ... \n", iR );
			#endif

			// Get distribution of next range
			rangeBlockIter->getDistribution( rangeBinData, &nb );

			// Find best match for the range distribution
			optRot = 0;
			optRef = false;
			bestMatch = matchRangeToDomains2( rangeBinData, &matchError, &optRot, &optRef, &muError, &sdError );

			#ifdef DEBUG_MODE
			fprintf( stderr, "%d-th range maps to %d-th domain with %g <%d, %d> error ...\n",
					 iR, bestMatch, matchError, optRot, (int)optRef );
			#endif

			// Add range information to domain
			#ifdef DEBUG_MODE
			fprintf( stderr, "Adding range information to the matched domain ... \n" );
			#endif
			typename list<FEL_domain<T> >::iterator tmpIter = domainList.begin();
			for( int iD=0; iD<bestMatch; iD++ )
				tmpIter++;
			tmpIter->addMatchingRange( bestMatch, matchError, optRot, optRef );

			// Add domain information to range
			#ifdef DEBUG_MODE
			fprintf( stderr, "Adding domain information to the matched range ... \n" );
			#endif
			rangeBlockIter->setMatchedDomainId( bestMatch );
			rangeBlockIter->setMatchError( matchError  );
			rangeBlockIter->setMatchTransformation( optRot, optRef );

			rangeBlockIter++;

		}// end for : scan through ranges

		#ifdef DEBUG_MODE
		fprintf( stderr, "Number of ranges: %d\n", rangeBlockList.size() );
		#endif

		// Clear up
		delete [] rangeBinData;
	}

	void
	encode_Regular_Distribution_WithBoundedError( double** rangeHistList, int nRange,
								 	 	 	 	  double errorBound )
	{
		int nb;
		double* rangeBinData = new double[nBin];
		int bestMatch;
		int optRot;
		bool optRef;
		double matchError, muError, sdError;
		float lowSub[3] = {0,0,0};
		float highSub[3] = {3, 3, 3};
		FEL_domain<T> newDomain( nDim, lowSub, highSub, 0, nBin );

		//#ifdef DEBUG_MODE
		fprintf( stderr, "Number of domain and range histograms: %d, %d\n", nTotalDomain, nRange );
		//#endif

		// Put all ranges in a list
		for( int iR=0; iR<nRange; iR++ )
		{
			// Copy range distribution
			memcpy( rangeBinData, rangeHistList[iR], sizeof(double)*nBin );

			// Create a range block and
			// set properties of the range block including distribution
			FEL_range<T> nextRange( 3, lowSub, highSub, nBin );
			nextRange.setDistribution( rangeBinData, nBin );

			// Add new range to the list
			rangeBlockList.push_back( nextRange );
		}

		//#ifdef DEBUG_MODE
		fprintf( stderr, "All ranges are in a list\n" );
		//#endif

		// Loop for encoding
		typename list< FEL_range<T> >::iterator rangeBlockIter = rangeBlockList.begin();
		for( int iR=0; iR<nRange; iR++ )
		{
			#ifdef DEBUG_MODE
			fprintf( stderr, "Finding best match for the %d-th range ... \n", iR );
			#endif

			// Get distribution of next range
			rangeBlockIter->getDistribution( rangeBinData, &nb );

			// Find best match for the range distribution
			optRot = 0;
			optRef = false;
			bestMatch = matchRangeToDomains2( rangeBinData, &matchError, &optRot, &optRef, &muError, &sdError );

			#ifdef DEBUG_MODE
			fprintf( stderr, "%d-th range maps to %d-th domain with %g <%d, %d> error ...\n",
					 iR, bestMatch, matchError, optRot, (int)optRef );
			#endif

			// Special case: if there was no match,
			// add the current range to the domain list
			if( matchError > errorBound )
			{
				//fprintf( stderr, "New domain needs to be added\n" );
				newDomain.setDistribution( rangeBinData, nBin );
				//fprintf( stderr, "here\n" );
				domainList.push_back( newDomain );
				//fprintf( stderr, "here2\n" );
				nTotalDomain = domainList.size();
				bestMatch = nTotalDomain - 1;
				optRot = 0;
				optRef = false;
				//fprintf( stderr, "here3\n" );
			}

			// Add range information to domain
			#ifdef DEBUG_MODE
			fprintf( stderr, "Adding range information to the matched domain ... \n" );
			#endif
			typename list<FEL_domain<T> >::iterator tmpIter = domainList.begin();
			for( int iD=0; iD<bestMatch; iD++ )
				tmpIter++;
			tmpIter->addMatchingRange( bestMatch, matchError, optRot, optRef );

			// Add domain information to range
			#ifdef DEBUG_MODE
			fprintf( stderr, "Adding domain information to the matched range ... \n" );
			#endif
			rangeBlockIter->setMatchedDomainId( bestMatch );
			rangeBlockIter->setMatchError( matchError  );
			rangeBlockIter->setMatchTransformation( optRot, optRef );

			rangeBlockIter++;

		}// end for : scan through ranges

		//#ifdef DEBUG_MODE
		fprintf( stderr, "Number of domains: %d\n", nTotalDomain );
		//#endif

		// Clear up
		delete [] rangeBinData;
	}

	void
	encode_Regular_Distribution_WithBoundedError2( double** rangeHistList, int nRange,
												   FEL_encodedhistogram* encodedHistList,
								 	 	 	 	   double errorThreshold,
								 	 	 	 	   int* domainUtilizationList,
								 	 	 	 	   int* nTotalHighErrorBin )
	{
		double* rangeFreqList = new double[nBin];
		int nb = 0, bestMatch = 0, optRot = 0, nHighErrorBin = 0;
		bool optRef = false;
		double matchError = 0;
		nUnMatchedRange = 0;
		(*nTotalHighErrorBin) = 0;
		int* highErrorBinIds = 0;
		double* highErrors = 0;


		#ifdef DEBUG_MODE
		fprintf( stderr, "Number of domain and range histograms: %d, %d\n", nTotalDomain, nRange );
		#endif

		// Loop for encoding
		for( int iR=0; iR<nRange; iR++ )
		{
			#ifdef DEBUG_MODE
			fprintf( stderr, "Finding best match for the %d-th range ... \n", iR );
			#endif

			// Get distribution of next range
			memcpy( rangeFreqList, rangeHistList[iR], sizeof(double)*nBin );

			// Find best match for the range distribution
			optRot = 0;
			optRef = false;
			bestMatch = matchRangeToDomains2( rangeFreqList, &matchError, &optRot, &optRef );

			// Keep track of utilized domains
			assert( bestMatch >= 0 );
			domainUtilizationList[bestMatch] ++;

			#ifdef DEBUG_MODE
			fprintf( stderr, "%d-th range maps to %d-th domain with %g <%d, %d> error ...\n",
					 iR, bestMatch, matchError, optRot, (int)optRef );
			#endif

			// Store mapping information to encoded histogram
			encodedHistList[iR].setMatchingDomainId( bestMatch );
			encodedHistList[iR].setRotation( optRot );
			encodedHistList[iR].setReflection( optRef );

			// Store high error bins
			if( errorThreshold != 0 && matchError > errorThreshold )
			{
				double* domainDist = new double[nBin];
				double* transformedDomain = new double[nBin];
				int* highErrorBinIds = new int[nBin];
				double* highErrors = new double[nBin];

				// Get domain distribution
				typename list<FEL_domain<T> >::iterator iter = domainList.begin();
				for( int i=0; i< bestMatch; i++ )
					iter++;
				iter->getDistribution( domainDist, &nb );

				// Transform domain
				memcpy( transformedDomain, domainDist, sizeof(double)*nBin );
				transformDistribution( transformedDomain, nBin, optRot, optRef );

				// Find high error bins and store errors
				int iP = 0;
				nHighErrorBin = 0;
				for( int i=0; i<nBin; i++ )
				{
					if( abs( transformedDomain[i] - rangeFreqList[i] ) > errorThreshold )
					{
						highErrorBinIds[iP] = i;
						highErrors[iP] = ( rangeFreqList[i] - transformedDomain[i] );

						iP++;
						nHighErrorBin ++;
						(*nTotalHighErrorBin) ++;
					}
				}

				// Store high error bins
				encodedHistList[iR].setNumHighErrorBin( nHighErrorBin );
				encodedHistList[iR].setHighErrorBinInfo( highErrorBinIds, highErrors );

				delete [] highErrorBinIds;
				delete [] highErrors;

			}// end if

		}// end for : scan through ranges

		#ifdef DEBUG_MODE
		fprintf( stderr, "Total number of unmatched range histograms: %d\n", nUnMatchedRange );
		fprintf( stderr, "Total number of high error bins: %d\n", (*nTotalHighErrorBin) );
		int nUsedDomain = 0;
		for( int i=0; i<nTotalDomain; i++ )
		{
			fprintf( stderr, "%d, ", domainUtilizationList[i] );
			if( domainUtilizationList[i] > 0 )
				nUsedDomain ++;
		}
		fprintf( stderr, "\n" );
		fprintf( stderr, "Total number of used domains: %d\n", nUsedDomain );
		#endif

		// Clear up
		delete [] rangeFreqList;
	}

	int
	storeHighErrorBins( double* rangeBinData, double* bestMatchDomainBinData,
						double errorBound,
						int optRot, bool optRef )
	{
		int nHighErrorBin = 0;

		for( int i=0; i<nBin; i++ )
		{
			if( fabs( rangeBinData[i] - bestMatchDomainBinData[i] ) > errorBound )
			{
				nHighErrorBin ++;
			}
		}

		return nHighErrorBin;
	}


	int
	matchRangeToDomains2( double* rangeFreqList, double *error, int* rotateAmount, bool* isReflected,
						  double* muerror = NULL, double* sderror = NULL )
	{
		int domainId = 0, bestMatchedDomainId = -1;
		double e = 0;
		double minError = 100;
		double domainFreqList[nBin];
		double bestDomainFreqList[nBin];
		int rotAmount = 0;
		bool isRef = false;

		for( typename list< FEL_domain<T> >::iterator domainIter = domainList.begin();
			 domainIter != domainList.end();
			 domainIter++ )
		{
			domainIter->getDistribution( domainFreqList, &nBin );
			#ifdef DEBUG_MODE
			double m = ITL_util<double>::Min( domainFreqList, nBin );
			double M = ITL_util<double>::Max( domainFreqList, nBin );
			fprintf( stderr, "%d-th Domain distribution min/max: %g, %g\n", domainID, m, M );
			#endif

			//e =  matchRangeToDomain2( rangeFreqList, iD, optTransformationMatrix );
			e =  FEL_util<float>::matchRangeToDomainByDistribution( distanceMeasureType,
																	rangeFreqList, domainFreqList,
																	nBin,
																	&rotAmount, &isRef );

			// Initialize min error
			if( domainId == 0 )	minError = e;

			if( e < minError )
			{
				minError = e;
				(*rotateAmount) = rotAmount;
				(*isReflected) = isRef;
				bestMatchedDomainId = domainId;
				memcpy( bestDomainFreqList, domainFreqList, sizeof( double )*nBin );
			}

			domainId ++;
		}

		// Optional analysis step: compute mean diffference and sd difference
		if( muerror != NULL )
			(*muerror) = FEL_util<double>::computeHistogramMeanError( rangeFreqList, bestDomainFreqList, nBin );
		if( sderror != NULL )
			(*sderror) = FEL_util<double>::computeHistogramSDError( rangeFreqList, bestDomainFreqList, nBin );


		*error = minError;
		return bestMatchedDomainId;

	}

	//list<ITL_distribution>*
	void
	identifySalientDistributions( float transformationSimilarityThreshold )
	{
		bool isRef;
		int nb, nRot;
		int nMatchedRange = 0;
		double domainFreqList[nBin];

		salientRangeDistributionList.clear();

		int domainID = 0;
		for( typename list<FEL_domain<T> >::iterator domIter = domainList.begin();
			 domIter != domainList.end();
			 domIter++ )

		{
			// Return and report matching ranges
			nMatchedRange = domIter->getNumMatchingRange();
			fprintf( stderr, "%d-th domain has %d matching ranges\n", domainID, nMatchedRange );
			if( nMatchedRange == 0 )
			{
				fprintf( stderr, "Continuing to next domain ...\n" );
				continue;
			}

			// Get domain distribution
			fprintf( stderr, "Scanning through matching range blocks ...\n" );
			domIter->getDistribution( domainFreqList, &nb );
			cout << "got domain distribution" << endl;

			int lastRot = -100;
			bool lastRef = false;
			int curRot = -100;
			double lastError, curError;
			bool curRef = false;
			list<int> tmpTrList;
			bool tmpRef;
			double avgDist[nBin];
			double trDist[nBin];

			memset( avgDist, 0, sizeof(double)*nBin );

			for( int i=0; i<nMatchedRange; i++ )
			{
				fprintf( stderr, "Continuing to next matching range ...\n" );
				int rID = domIter->getMatchingRangeInformation(i, &curError, &curRot, &curRef );

				if( i == 0 )
				{
					lastError = curError;
					lastRot = curRot;
					lastRef = curRef;
					tmpTrList.push_back( curRot );
					tmpRef = curRef;
					continue;
				}
				else
				{
					if( (curError - lastError) <= transformationSimilarityThreshold )
					{
						// Transform domain distribution to approximate range distribution
						memcpy( trDist, domainFreqList, sizeof(double)*nBin );
						transformDistribution( trDist, nBin, curRot, curRef );

						// Update representative distribution
						ITL_util<double>::addArrays( avgDist, trDist, avgDist, nBin );

						// Normalize representative distribution
						float sum = 0;
						for( int i=0; i<nBin; i++ )
							sum += avgDist[i];
						for( int i=0; i<nBin; i++ )
							avgDist[i] = avgDist[i] / sum;
					}
					else
					{
						// Store and save the current stack
						ITL_distribution salientDist( nBin );
						salientDist.setFrequencies( avgDist );
						salientRangeDistributionList.push_back( salientDist );

						// Transform domain distribution to match range
						memcpy( trDist, domainFreqList, sizeof(float)*nBin );
						transformDistribution( trDist, nBin, curRot, curRef );

						memcpy( avgDist, trDist, sizeof(double)*nBin );

					}// end else : curerror - minerror

					lastError = curError;
					lastRot = curRot;
					lastRef = curRef;

				}// end outer else

			}// end inner for

			domainID ++;

		}// end outer for

		fprintf( stderr, "%d salient distributions obtained from the data ... \n", salientRangeDistributionList.size() );

	}// End function


	list<FEL_range<T> >
	searchDistribution( ITL_distribution* queryDistribution )
	{
		list<FEL_range<T> > returnList;
		double qFreqList[nBin];
		float qEntropy;
		float qRange[2];
		float delH;
		list<int> searchResultDomainList;
		int rangeBounds[6];
		float lowF[3], highF[3];
		int nReturnedRange;
		FEL_range<T>* nextRange = NULL;

		double returnedRangeFreqList[nBin];
		double matchError;
		int nMatchedRange;
		//int* rangeIdList = NULL;
		int rangeIdList[2000];
		int nb;

		// Compute entropy of the query distribution
		queryDistribution->getFrequencies( qFreqList );
		qEntropy = ITL_entropycore::computeEntropy_HistogramBased2( qFreqList, nBin, false );
		fprintf( stderr, "Entropy of query distribution: %g\n", qEntropy );

		// Store distribution for Matlab
		//for( int iB = 0; iB<nBin; iB++ )
		//	queryDist[iB] = (double)qFreqList[iB];

		// Create a query entropy range
		delH = ( 0.05f * ( log(nBin) / log(2) ) );
		qRange[0] = qEntropy - delH;
		qRange[1] = qEntropy + delH;
		fprintf( stderr, "Entropy range for querying: <%g, %g>\n", qRange[0], qRange[1] );

		// Step 1: Return domains within query range
		searchResultDomainList = domainTree.searchTree( domainTree.getRoot(), qRange[0], qRange[1] );
		#ifdef DEBUG_MODE
		fprintf( stderr, "Domain search result:\n" );
		for( list<int>::iterator it = searchResultDomainList.begin();
			 it != searchResultDomainList.end();
			 it++ )
		{
			fprintf( stderr, "%d\t", (*it) );
		}
		fprintf( stderr, "\n" );
		#endif
		cout << "c" << endl;

		// For each returned domain, find the corresponding ranges
		nReturnedRange = 0;
		for( list<int>::iterator it = searchResultDomainList.begin();
			 it != searchResultDomainList.end();
			 it++ )
		{
			// Retrieve next domain
			int nextDomain = (*it);
			typename list<FEL_domain<T> >::iterator tmpIter = domainList.begin();
			for( int iD=0; iD<nextDomain; iD++ )
				tmpIter++;

			// Return and report matching ranges
			nMatchedRange = tmpIter->getNumMatchingRange();
			fprintf( stderr, "%d-th domain has %d matching ranges\n", nextDomain, nMatchedRange );

			if( nMatchedRange == 0 )
				continue;
			else
			{
				cout << nMatchedRange << endl;
				assert( nMatchedRange <= 2000 );
				//rangeIdList = new int[nMatchedRange];
				tmpIter->getMatchingRange( rangeIdList );
				cout << nMatchedRange << " " << tmpIter->getNumMatchingRange() << endl;

				// Scan though the range blocks
				// that match with the returned domain
				for( int iR=0; iR<nMatchedRange; iR++ )
				{
					if( encodingType == REGULAR )
					{
						//cout << "here" << endl;
						nextRange = searchRangeByID( rangeIdList[iR] );
						assert( nextRange != NULL );
						//cout << "here2" << endl;

					}
					else
					{
						nextRange = rangeBlockTree.searchTreeByRangeID( rangeBlockTree.getRoot(), rangeIdList[iR] );
					}

					// Copy the returned range information
					// and add it to the return list
					if( nextRange != NULL )
					{
						nextRange->getBounds( rangeBounds );
						printf( "%d-th returned range: <%d %d %d> <%d %d %d>\n", rangeIdList[iR],
												rangeBounds[0], rangeBounds[1], rangeBounds[2],
												rangeBounds[3], rangeBounds[4], rangeBounds[5] );

						nextRange->getDistribution( returnedRangeFreqList, &nb  );

						// Compute error between query and
						// the original (untransformed) range block distribution
						matchError = FEL_util<double>::computeHistogramMatchError( distanceMeasureType,
																				   qFreqList, returnedRangeFreqList,
																				   nBin );
						printf( "Match error between query and original range: %g\n", matchError );

						// Store selected range information
						// to return range
						lowF[0] = rangeBounds[0]; lowF[1] = rangeBounds[1]; lowF[2] = rangeBounds[2];
						highF[0] = rangeBounds[3]; highF[1] = rangeBounds[4]; highF[2] = rangeBounds[5];
						FEL_range<T> returnRange( nDim, lowF, highF, nBin );
						returnRange.setDistribution( returnedRangeFreqList, nBin );
						returnRange.setMatchError( matchError );

						// Add return range to return list
						//returnList->push_back( returnRange );
						returnList.push_back( returnRange );

						nReturnedRange ++;

					}// end if

					//delete [] rangeIdList;

				}// end if-else

			}// end for : iR

		}// end for : domain iterator

		return returnList;

	}// End function

	void
	searchDistribution2( ITL_distribution* queryDistribution, list<FEL_range<T> >* returnList )
	{
		double qFreqList[nBin];
		float qEntropy;
		float qRange[2];
		float delH;
		list<int> searchResultDomainList;
		int rangeBounds[6];
		float lowF[3], highF[3];
		int nReturnedRange;
		FEL_range<T>* nextRange = NULL;
		FEL_range<T> returnRange;
		double returnedRangeFreqList[nBin];
		double matchError;
		int nMatchedRange;
		int* rangeIdList = NULL;
		int nb;

		cout << "a" << endl;

		// Compute entropy of the query distribution
		queryDistribution->getFrequencies( qFreqList );
		qEntropy = ITL_entropycore::computeEntropy_HistogramBased2( qFreqList, nBin, false );
		fprintf( stderr, "Entropy of query distribution: %g\n", qEntropy );

		cout << "b" << endl;

		// Store distribution for Matlab
		//for( int iB = 0; iB<nBin; iB++ )
		//	queryDist[iB] = (double)qFreqList[iB];

		// Create a query entropy range
		delH = ( 0.1f * ( log(nBin) / log(2) ) );
		qRange[0] = qEntropy - delH;
		qRange[1] = qEntropy + delH;
		fprintf( stderr, "Entropy range for querying: <%g, %g>\n", qRange[0], qRange[1] );

		// Step 1: Return domains within query range
		searchResultDomainList = domainTree.searchTree( domainTree.getRoot(), qRange[0], qRange[1] );
		#ifdef DEBUG_MODE
		fprintf( stderr, "Domain search result:\n" );
		for( list<int>::iterator it = searchResultDomainList.begin();
			 it != searchResultDomainList.end();
			 it++ )
		{
			fprintf( stderr, "%d\t", (*it) );
		}
		fprintf( stderr, "\n" );
		#endif

		// For each returned domain, find the corresponding ranges
		nReturnedRange = 0;

		for( list<int>::iterator it = searchResultDomainList.begin();
			 it != searchResultDomainList.end();
			 it++ )
		{
			// Retrieve next domain
			int nextDomain = (*it);
			typename list<FEL_domain<T> >::iterator tmpIter = domainList.begin();
			for( int iD=0; iD<nextDomain; iD++ )
				tmpIter++;

			// Return and report matching ranges
			nMatchedRange = tmpIter->getNumMatchingRange();
			fprintf( stderr, "%d-th domain has %d matching ranges\n", nextDomain, nMatchedRange );

			if( nMatchedRange == 0 )
				continue;
			else
			{
				cout << nMatchedRange << endl;
				rangeIdList = new int[nMatchedRange];
				tmpIter->getMatchingRange( rangeIdList );
				cout << nMatchedRange << " " << tmpIter->getNumMatchingRange() << endl;

				// Scan though the range blocks
				// that match with the returned domain
				for( int iR=0; iR<nMatchedRange; iR++ )
				{
					if( encodingType == REGULAR )
					{
						cout << "here" << endl;
						nextRange = searchRangeByID( rangeIdList[iR] );
						assert( nextRange != NULL );
						cout << "here2" << endl;

					}
					else
					{
						nextRange = rangeBlockTree.searchTreeByRangeID( rangeBlockTree.getRoot(), rangeIdList[iR] );
					}

					// Copy the returned range information
					// and add it to the return list
					if( nextRange != NULL )
					{
						nextRange->getBounds( rangeBounds );
						printf( "%d-th returned range: <%d %d %d> <%d %d %d>\n", rangeIdList[iR],
												rangeBounds[0], rangeBounds[1], rangeBounds[2],
												rangeBounds[3], rangeBounds[4], rangeBounds[5] );

						nextRange->getDistribution( returnedRangeFreqList, &nb  );

						// Compute error between query and
						// the original (untransformed) range block distribution
						matchError = computeMatchError2( qFreqList, returnedRangeFreqList  );
						printf( "Match error between query and original range: %g\n", matchError );

						// Store selected range information
						// to return range
						lowF[0] = rangeBounds[0]; lowF[1] = rangeBounds[1]; lowF[2] = rangeBounds[2];
						highF[0] = rangeBounds[3]; highF[1] = rangeBounds[4]; highF[2] = rangeBounds[5];
						returnRange.initialize( nDim, lowF, highF );
						returnRange.setDistribution( returnedRangeFreqList, nBin );
						returnRange.setMatchError( matchError );

						// Add return range to return list
						returnList->push_back( returnRange );

						nReturnedRange ++;

					}// end if

					delete [] rangeIdList;

				}// end if-else

			}// end for : iR

		}// end for : domain iterator

	}// End function

	FEL_range<T>*
	searchRangeByID( int rangeID )
	{
		typename list<FEL_range<T> >::iterator rangeIter = rangeBlockList.begin();
		for( int iR=0; iR<rangeID; iR++ )
			rangeIter ++;

		// Is there a better/safer way to return pointers?
		return &(*rangeIter);

	}// End function


	void
	estimateDistributionInARegion( int* queryBounds, double* estDist )
	{
		assert( estDist!= NULL );
		memset( estDist, 0, sizeof(double)*nBin );
		int queryBoundsLocal[6];
		double estDistLocal[nBin];
		memset( estDistLocal, 0, sizeof(double)*nBin );

		// Return all zero if query region has size zero or negative
		if( queryBounds[0] > queryBounds[3] ||
			queryBounds[1] > queryBounds[4] ||
			queryBounds[2] > queryBounds[5] )
			return;

		if( encodingType == REGULAR )
		{
			memcpy( queryBoundsLocal, queryBounds, sizeof(int)*6 );
			estimateDistributionInARegion_Regular( queryBoundsLocal, estDistLocal );
			memcpy( estDist, estDistLocal, sizeof(double)*nBin );

		} // end if : regular encoding
		else
		{
			memcpy( queryBoundsLocal, queryBounds, sizeof(int)*6 );
			estimateDistributionInARegion_Treebased( queryBoundsLocal, estDistLocal );
			#ifdef DEBUG_MODE
			double m = ITL_util<double>::Min( estDist, nBin );
			double M = ITL_util<double>::Max( estDist, nBin );
			cout << "m1: " << m << " " << "M1: " << M << endl;
			#endif
			memcpy( estDist, estDistLocal, sizeof(double)*nBin );

		}// end else : tree based encoding

		#ifdef USING_PTHREAD
		pthread_exit(NULL);
		#endif

	}// End function

	void
	estimateDistributionInARegion_Regular( int* queryBounds, double* estDist )
	{
		int nb;
		int matchedDomainID;
		int nRot;
		bool isRef;
		double matchedDomainFreqList[nBin];
		int rangeBounds[6];

		// Traverse the list of ranges and return
		// any range that overlaps with the query
		for( typename list<FEL_range<T> >::iterator rangeIter = rangeBlockList.begin();
			 rangeIter != rangeBlockList.end();
			 ++rangeIter )
		{
			// Get range bounds
			rangeIter->getBounds( rangeBounds );

			// If range overlaps the query,
			// Get distribution of matching domain
			// and add to the estimated distribution
			if( FEL_util<int>::isOverlapping( queryBounds, rangeBounds ) )
			{
				// Get ID of matching domain
				matchedDomainID = rangeIter->getMatchedDomainId();

				//#ifdef DUBUG_MODE
				fprintf( stderr, "Query [<%d,%d,%d>-<%d,%d,%d>] overlaps with range [<%d,%d,%d>-<%d,%d,%d>]\n",
						 queryBounds[0], queryBounds[1], queryBounds[2],
						 queryBounds[3], queryBounds[4], queryBounds[5],
						 rangeBounds[0], rangeBounds[1], rangeBounds[2],
						 rangeBounds[3], rangeBounds[4], rangeBounds[5] );
				fprintf( stderr, "Matched Domain: %d\n", matchedDomainID );
				//#endif

				// Get distribution and transformation of matching domain
				typename list<FEL_domain<T> >::iterator domIter = domainList.begin();
				for( int iD=0; iD<matchedDomainID; iD++ )
					domIter++;
				domIter->getDistribution( matchedDomainFreqList, &nb );
				rangeIter->getMatchTransformation( nRot, isRef );

				// Transform domain distribution to match range
				transformDistribution( matchedDomainFreqList, nBin, nRot, isRef );

				// Update estimated distribution
				ITL_util<double>::addArrays( estDist, matchedDomainFreqList, estDist, nBin );

			}
		}

		// Normalize estimated distribution
		double sum = 0;
		for( int i=0; i<nBin; i++ )
			sum += estDist[i];
		for( int i=0; i<nBin; i++ )
			estDist[i] = estDist[i] / sum;

	}// End function

	void
	estimateDistributionInARegion_Treebased( int* queryBounds, double* estDist )
	{
		int nb;
		int matchedDomainID;
		int nRot;
		bool isRef;
		double matchedDomainFreqList[nBin];
		int rangeBounds[6];
		float queryBoundsLocal[6];
		list<FEL_range<T> > overlappingLeafRangeList;
		//FEL_range<T>* leafRange = NULL;
		double estDistLocal[nBin];
		memset( estDistLocal, 0, sizeof(double)*nBin );

		// Get a list of leaf nodes (range blocks)
		// that overlaps with the query
		for( int i=0; i<6; i++ )
			queryBoundsLocal[i] = queryBounds[i];
		overlappingLeafRangeList = rangeBlockTree.searchLeafNodesByRegion( rangeBlockTree.getRoot(), queryBoundsLocal );
		#ifdef DUBUG_MODE
		fprintf( stderr, "%d leaf nodes overlap with the query ...\n", overlappingLeafRangeList.size() );
		#endif

		// Scan through the list
		// of returned ranges and collect distributions
		//for( typename list<FEL_range<T>* >::iterator rangeIter = overlappingLeafRangeList.begin();
		double qRatio = 1;
		double totalWeight = 0;
		for( typename list<FEL_range<T> >::iterator rangeIter = overlappingLeafRangeList.begin();
			 rangeIter != overlappingLeafRangeList.end();
			 ++rangeIter )
		{
			// Get next leaf range
			//leafRange = *rangeIter;

			// Get range bounds
			rangeIter->getBounds( rangeBounds );
			//leafRange->getBounds( rangeBounds );

			// Compute query portion size / range block size
			//qRatio = (double)FEL_util<int>::getOverlapSize( queryBounds, rangeBounds ) / (double)rangeIter->getSize();
			//fprintf( stderr, "Q-Ratio: %g\n", qRatio );

			#ifdef DUBUG_MODE
			fprintf( stderr, "Query [<%d,%d,%d>-<%d,%d,%d>] overlaps with range [<%d,%d,%d>-<%d,%d,%d>]\n",
					 queryBounds[0], queryBounds[1], queryBounds[2],
					 queryBounds[3], queryBounds[4], queryBounds[5],
					 rangeBounds[0], rangeBounds[1], rangeBounds[2],
					 rangeBounds[3], rangeBounds[4], rangeBounds[5] );
			fprintf( stderr, "Matched Domain: %d\n", matchedDomainID );
			#endif

			// Get ID of matching domain
			// If their is no good matching domain,
			// then obtain the
			if( rangeIter->	getMatQualityFlag() )
			{
				matchedDomainID = rangeIter->getMatchedDomainId();
				//matchedDomainID = leafRange->getMatchedDomainId();

				// Get distribution and transformation of matching domain
				typename list<FEL_domain<T> >::iterator domIter = domainList.begin();
				for( int iD=0; iD<matchedDomainID; iD++ )
					domIter++;
				domIter->getDistribution( matchedDomainFreqList, &nb );
				rangeIter->getMatchTransformation( nRot, isRef );
				//leafRange->getMatchTransformation( nRot, isRef );en
				//cout << nRot << " " << isRef << endl;

				// Transform domain distribution to match range
				transformDistribution( matchedDomainFreqList, nBin, nRot, isRef );
			}
			else
				rangeIter->getDistribution( matchedDomainFreqList, &nb );

			// Update estimated distribution
			//if( qRatio < 1 )
			//{
			//	FEL_util<double>::scaleDistribution( matchedDomainFreqList, nBin, qRatio );
			//	totalWeight += qRatio;
			//}
			//else
			//	totalWeight += 1;

			ITL_util<double>::addArrays( estDistLocal, matchedDomainFreqList, estDistLocal, nBin );
			#ifdef DEBUG_MODE
			double m = ITL_util<double>::Min( estDistLocal, nBin );
			double M = ITL_util<double>::Max( estDistLocal, nBin );
			cout << "m1: " << m << " " << "M1: " << M << endl;
			#endif

		}// end for : overlappingLeafRangeList

		// Normalize estimated distribution
		double sum = 0;
		for( int i=0; i<nBin; i++ )
			sum += estDistLocal[i];
		for( int i=0; i<nBin; i++ )
			estDistLocal[i] = estDistLocal[i] / sum; //(sum*totalWeight);

		memcpy( estDist, estDistLocal, sizeof(double)*nBin );

	}// End function

	void
	decode_Regular_AllDistributions( int nRange,
									 FEL_encodedhistogram* encodedHistList,
									 double** decodedDist )
	{
		int nb = 0, matchingDomainId = 0, optRot = 0, nHighErrorBin = 0;
		bool optRef;
		long nDataElement;
		double* matchedDomainFreqList = new double[nBin];
		double* freqList = new double[nBin];
		int rangeBounds[6];
		double sum = 0;

		// Traverse the list of ranges and return
		// any range that overlaps with the query
		for( int iR=0; iR<nRange; iR++ )
		{
			// Get information from encoder
			matchingDomainId = encodedHistList[iR].getMatchingDomainId();
			nDataElement = encodedHistList[iR].getNumDataElement();
			nHighErrorBin = encodedHistList[iR].getNumHighErrorBin();
			optRot = encodedHistList[iR].getRotation();
			optRef = encodedHistList[iR].getReflection();

			// Get domain distribution with specified ID
			typename list<FEL_domain<T> >::iterator domainIter = domainList.begin();
			for( int iD=0; iD<matchingDomainId; iD++ )
				domainIter++;
			domainIter->getDistribution( freqList, &nb );

			#ifdef DEBUG_MODE
			fprintf( stderr, "In decoder: matching domain id: %d\n", matchingDomainId );
			fprintf( stderr, "In decoder: number of elements: %d\n", nDataElement );
			fprintf( stderr, "In decoder: rotation and reflection: %d %d\n", optRot, optRef );
			fprintf( stderr, "In decoder: Number of high error bins: %d\n", nHighErrorBin );
			fprintf( stderr, "In decoder: matching domain distribution:\n" );
			for( int i = 0; i<nBin; i++ )
				fprintf( stderr, "%g, ", freqList[i] );
			fprintf( stderr, "\n" );
			#endif

			// Transform domain to match range
			transformDistribution( freqList, nBin, optRot, optRef );
			#ifdef DEBUG_MODE
			fprintf( stderr, "In decoder: transformed domain distribution:\n" );
			for( int i = 0; i<nBin; i++ )
				fprintf( stderr, "%g, ", freqList[i] );
			fprintf( stderr, "\n" );
			#endif

			// Further refine the high error bins
			if( nHighErrorBin > 0 )
			{
				int* highErrorBinIds = new int[nHighErrorBin];
				double* highErrors = new double[nHighErrorBin];

				// Get stored error information
				encodedHistList[iR].getHighErrorBinInfo( highErrorBinIds, highErrors );
				#ifdef DEBUG_MODE
				fprintf( stderr, "In decoder: Error Info:\n" );
				for( int i = 0; i<nHighErrorBin; i++ )
					fprintf( stderr, "%d: %g\n", highErrorBinIds[i], highErrors[i] );
				#endif

				// Use information to improve decoding
				correctHighErrorBins( freqList, nHighErrorBin, highErrorBinIds, highErrors );
				#ifdef DEBUG_MODE
				fprintf( stderr, "In decoder: error corrected decoded distribution:\n" );
				for( int i = 0; i<nBin; i++ )
					fprintf( stderr, "%g, ", freqList[i] );
				fprintf( stderr, "\n" );
				#endif

				// Renormalize
				double sumFreq = ITL_util<double>::sum( matchedDomainFreqList, nBin );
				ITL_util<double>::divideArrayScalar( matchedDomainFreqList, sumFreq, nBin );

				// Free
				delete [] highErrorBinIds;
				delete [] highErrors;
			}

			// Multiply number of elements with normalized distribution
			if( nDataElement != 0 )
			{
				for( int i=0; i<nBin; i++ )
					freqList[i] = freqList[i] * nDataElement;
			}

			// Copy decoded distribution
			memcpy( decodedDist[iR], freqList, sizeof(double)*nBin );
		}// end for : scan through ranges


		// Clear up
		delete [] matchedDomainFreqList;
		delete [] freqList;


	}// end function

	void
	decode_Distribution( FEL_encodedhistogram* cbe, double* decodedFreqList )
	{
		int matchedDomainID, nb, nRot, nHighError;
		bool isRef;
		double* matchedDomainFreqList = new double[nBin];

		// Get required information from encoded histogram
		matchedDomainID = cbe->getMatchingDomainId();
		nRot = cbe->getRotation();
		isRef = cbe->getReflection();
		nHighError = cbe->getNumHighErrorBin();

		int* binIds = new int[nHighError];
		double* highErrors = new double[nHighError];
		cbe->getHighErrorBinInfo( binIds, highErrors );

		#ifdef DUBUG_MODE
		fprintf( stderr, "Matched Domain: %d\n", matchedDomainID );
		#endif

		// Get matching domain distribution
		typename list<FEL_domain<T> >::iterator domIter = domainList.begin();
		for( int iD=0; iD<matchedDomainID; iD++ )
			domIter++;
		domIter->getDistribution( matchedDomainFreqList, &nb );

		// Transform domain distribution to match range
		transformDistribution( matchedDomainFreqList, nBin, nRot, isRef );

		// Correct high error bins
		correctHighErrorBins( matchedDomainFreqList, nHighError, binIds, highErrors );

		// Renormalize
		double sumFreq = ITL_util<double>::sum( matchedDomainFreqList, nBin );
		ITL_util<double>::divideArrayScalar( matchedDomainFreqList, sumFreq, nBin );

		// Copy decoded histogram
		memcpy( decodedFreqList, matchedDomainFreqList, sizeof(double)*nBin );

		// Clear up
		delete [] matchedDomainFreqList;
		delete [] binIds;
		delete [] highErrors;



	}// end function

	void
	transformDistribution( double* freqList, int nBin, int nRot, bool isRef )
	{
		if( nRot == 0 && isRef == false )
			return;

		// Create buffer
		double *freqListTemp = new double[nBin];

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
			/*
			double buffer;
			double buffer2[nBin-1];
			for( int iT = 0; iT<nRot; iT++ )
			{
				// Place rightmost to buffer
				buffer = freqListTemp[nBin-1];

				// Move everybody else by 1 bit
				//memcpy( buffer2, freqListTemp, sizeof(float)*(nBin-1) );
				//memcpy( freqListTemp+1, buffer2, sizeof(float)*(nBin-1) );
				for( int j=nBin-1; j>=1; j-- )
					freqListTemp[j] = freqListTemp[j-1];

				// Fill back the leftmost
				freqListTemp[0] = buffer;
			} // end for : iT
			*/
			double buffer[nRot];
			double buffer2[nBin-nRot];
			memcpy( buffer, freqListTemp, sizeof(double)*nRot );
			memcpy( buffer2, freqListTemp+nRot, sizeof(double)*(nBin-nRot) );
			memcpy( freqListTemp, buffer2, sizeof(double)*(nBin-nRot) );
			memcpy( freqListTemp+(nBin-nRot), buffer, sizeof(double)*nRot );

		}

		memcpy( freqList, freqListTemp, sizeof(double)*nBin );

		// Free temporary resources
		delete [] freqListTemp;
	}

	void
	correctHighErrorBins( double* freqList, int nHighErrorBin, int* highErrorBinId, double* highError )
	{
		if( nHighErrorBin == 0 )
			return;

		int id = -1;
		for( int i=0; i<nHighErrorBin; i++ )
		{
			id = highErrorBinId[i];
			freqList[id] = freqList[id] + highError[i];
		}
	}// end function

	void
	filterDomainBlockExtents( int lod, int* nReturnedDomains, float* domainExtentList )
	{
		assert( domainExtentList != NULL );
		(*nReturnedDomains) = 0;

		int bounds[6];
		int iD = 0;
		for( typename list<FEL_domain<T> >::iterator tmpIter = domainList.begin();
			 tmpIter != domainList.end();
			 tmpIter++ )
		{
			if( tmpIter->getDataLoD() != lod )
				continue;

			tmpIter->getBounds( bounds );
			for( int i=0; i<6; i++ )
				domainExtentList[iD*6+i] = (float)bounds[i];

			(*nReturnedDomains) ++;
			iD++;

		}// end for
	}

	//FEL_range<VECTOR3> getRangeBlock( int id );
	void
	getRangeBlockDistribution( int rangeid, double* fList, int* matchingDomainID,
								double* error, int* optRot, bool* optRef )
	{
		assert( fList != NULL );
		double tmp[nBin];
		int nb;

		if( encodingType == REGULAR )
		{
			fprintf( stderr, "Searching range block list for the range block with ID: %d\n", rangeid );
			// Linearly scan the range list to find the queried one
			// if regular encoding was used
			typename list<FEL_range<T> >::iterator tmpIter = rangeBlockList.begin();
			for( int iD=0; iD<rangeid; iD++ )
				tmpIter++;

			// Get distribitution
			tmpIter->getDistribution( tmp, &nb );

			// Get matching domain ID
			(*matchingDomainID) = tmpIter->getMatchedDomainId();

			// Get Match error
			(*error) = tmpIter->getMatchError();

			// Get optimal transformation
			tmpIter->getMatchTransformation( *optRot, *optRef );

			// Copy
			memcpy( fList, tmp, sizeof(float)*nBin );
		}
		else
		{
			fprintf( stderr, "Searching range block tree for the range block with ID: %d\n", rangeid );
			// Search quadtree/octree to find the queried one
			// if hierarchical encoding was used
			FEL_range<T>* rangePtr = rangeBlockTree.searchTreeByRangeID( rangeBlockTree.getRoot(), rangeid );

			//cout << "1" << endl;
			assert( rangePtr != NULL );

			// Get distribitution
			rangePtr->getDistribution( tmp, &nb );

			//cout << "2" << endl;

			// Get matching domain
			(*matchingDomainID) = rangePtr->getMatchedDomainId();

			// Get Match error
			(*error) = rangePtr->getMatchError();

			// Get optimal transformation
			rangePtr->getMatchTransformation( *optRot, *optRef );

			//cout << "3" << endl;

			// Copy
			memcpy( fList, tmp, sizeof(double)*nBin );

		}

		//double m = ITL_util<double>::Min( tmp, nBin );
		//double M = ITL_util<double>::Max( tmp, nBin );
		//cout << "m1: " << m << " " << "M1: " << M << endl;

	}// end function

	void filterRangeBlockExtents( int lod, int* nReturnedRanges, float* rangeExtentList )
	{
		assert( rangeExtentList != NULL );
		list< FEL_range<T>* > rangeAtLodList;
		int bounds[6];
		int iD = 0;

		if( encodingType == REGULAR )
		{
			fprintf( stderr, "All range blocks at the same level, so returning all ...\n" );
			(*nReturnedRanges) = nTotalRange;

			for( typename list<FEL_range<T> >::iterator tmpIter = rangeBlockList.begin();
				 tmpIter != rangeBlockList.end();
				 tmpIter++ )
			{
				tmpIter->getBounds( bounds );
				for( int i=0; i<6; i++ )
					rangeExtentList[iD*6+i] = (float)bounds[i];

				iD++;
			}// end for
		}
		else
		{
			fprintf( stderr, "Searching range block tree for the range blocka at level: %d\n", lod );
			if( lod >= 0 )
				rangeAtLodList = rangeBlockTree.searchTreeByLevel( rangeBlockTree.getRoot(), lod );
			else if( lod == -1 )
				rangeAtLodList = rangeBlockTree.getLeafNodes( rangeBlockTree.getRoot() );

			(*nReturnedRanges) = rangeAtLodList.size();
			fprintf( stderr, "Number of range block levels at this level: %d\n", (*nReturnedRanges) );
			if( (*nReturnedRanges) == 0 )
				return;

			fprintf( stderr, "Scanning through the returned lists and copying extents ...\n" );
			for( typename list<FEL_range<T>* >::iterator tmpIter = rangeAtLodList.begin();
				 tmpIter != rangeAtLodList.end();
				 tmpIter++ )
			{
				FEL_range<T>* tmp = *tmpIter;
				tmp->getBounds( bounds );
				for( int i=0; i<6; i++ )
					rangeExtentList[iD*6+i] = (float)bounds[i];

				iD++;
			}// end for


		}// end else
	}// end function

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Get-set functions
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	setHistogram( ITL_histogram *hist, int nbin, int mappingtype )
	{
		histogram = hist;
		nBin = nbin;
		histMappingType = mappingtype;
	}// end function

	FEL_spacetree<T>*
	getSpaceTree()
	{
		return ( &rangeBlockTree );
	}// end function

	int
	getNumTotalRange()
	{
		return nTotalRange;
	}// end function

	int
	getNumTotalDomain()
	{
		return nTotalDomain;
	}// end function

	int
	getDatatype()
	{
		return dataType;
	}// end function

	FEL_domain<T>
	getDomainBlock( int id )
	{
		typename list<FEL_domain<T> >::iterator tmpIter = domainList.begin();
		for( int iD=0; iD<id; iD++ )
			tmpIter++;

		double tmp[nBin];
		int nb;
		tmpIter->getDistribution( tmp, &nb );

		FEL_domain<T> retDomain = (*tmpIter);
		return retDomain;

	}// end function

	void
	getDomainBlockExtent( int id, float *extent )
	{
		assert( extent != NULL );
		int localExtent[6];

		typename list<FEL_domain<T> >::iterator tmpIter = domainList.begin();
		for( int iD=0; iD<id; iD++ )
			tmpIter++;

		tmpIter->getBounds( localExtent );
		if( domainCreationType == SINGLERES )
		{
			for( int i=0; i<6; i++ )
				extent[i] = localExtent[i]*2;
		}
		else
		{
			int lod = tmpIter->	getDataLoD();
			for( int i=0; i<6; i++ )
				extent[i] = localExtent[i]*pow( 2, lod+1 );
		}

		// Special case for 2D vector fields
		if( localExtent[5] == 1 )
		{
			extent[2] = localExtent[2];
			extent[5] = localExtent[5];
		}
	}// end function

	void
	getDomainBlockDistribution( int id, double *fList )
	{
		assert( fList != NULL );

		typename list<FEL_domain<T> >::iterator tmpIter = domainList.begin();
		for( int iD=0; iD<id; iD++ )
			tmpIter++;

		double tmp[nBin];
		int nb;
		tmpIter->getDistribution( tmp, &nb );
		memcpy( fList, tmp, sizeof(double)*nBin );

		#ifdef DEBUG_MODE
		double m = ITL_util<double>::Min( tmp, nBin );
		double M = ITL_util<double>::Max( tmp, nBin );
		cout << "m1: " << m << " " << "M1: " << M << endl;
		#endif

	}// end function

	void
	getRangeBlockExtent( int rangeid, float* extent )
	{
		int localExtent[6];
		assert( extent != NULL );

		if( encodingType == REGULAR )
		{
			fprintf( stderr, "Searching range block list for the range block with ID: %d\n", rangeid );
			// Linearly scan the range list to find the queried one
			// if regular encoding was used
			typename list<FEL_range<T> >::iterator tmpIter = rangeBlockList.begin();
			for( int iD=0; iD<rangeid; iD++ )
				tmpIter++;

			// Get bounds
			tmpIter->getBounds( localExtent );
			for( int i=0; i<6; i++ )
				extent[i] = localExtent[i];

		}
		else
		{
			fprintf( stderr, "Searching range block tree for the range block with ID: %d\n", rangeid );
			// Search quadtree/octree to find the queried one
			// if hierarchical encoding was used
			FEL_range<T>* rangePtr = rangeBlockTree.searchTreeByRangeID( rangeBlockTree.getRoot(), rangeid );

			// Get bounds
			rangePtr->getBounds( localExtent );
			for( int i=0; i<6; i++ )
				extent[i] = localExtent[i];
		}
	}// end function

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Printing, plotting and debugging related functions
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	toggleMatlabPlotting( bool flag )
	{
		if( flag )
		{
			// Do nothing if plotter is already running
			if( plotter != NULL )	return;
			else
			{
				plotter = new FEL_matlab();
				plotter->initEngine();

				matlab_plotting_on = true;
			}
		}
		else
		{
			// Do nothing if plotter is already off
			if( plotter == NULL )	return;
			else
			{
				plotter->closeEngine();
				delete [] plotter;

				cout << "114" << endl;
				matlab_plotting_on = false;
			}
		}
	}// end function

	void
	plotEncodingResutls()
	{
		// Draw Matlab plots
		cout << "11" << endl;
		//if( matlab_plotting_on )
		//{
			cout << "12" << endl;
			assert( plotter != NULL );
			plotter->createMatlabVariables();
			plotter->displayDomainUtilization_Matlab();
			plotter->displayMatchingPairEntropyLevel_Matlab();
			plotter->displayMatchErrorDistribution_Matlab();
			plotter->displayUnMatchedRanges_Matlab();

		//}
	}// end function

	void
	printDomains()
	{
		int domainBounds[6];
		int domainID = 0;

		for( typename list< FEL_domain<T> >::iterator domainIter = domainList.begin();
			 domainIter != domainList.end();
			 domainIter++ )
		{
			domainIter->getBounds( domainBounds );
			printf( "%d-th domain: <%d %d %d> <%d %d %d>, Entropy: %g\n", domainID,
									domainBounds[0], domainBounds[1], domainBounds[2],
									domainBounds[3], domainBounds[4], domainBounds[5],
									domainIter->getEntropy() );
			domainID ++;
		}

	}// end function

	void
	printRanges()
	{
		int rangeBounds[6];
		int rangeId = 0;

		for( typename list< FEL_range<T> >::iterator rangeBlockIter = rangeBlockList.begin();
			 rangeBlockIter != rangeBlockList.end();
			 rangeBlockIter++ )
		{
			rangeBlockIter->getBounds( rangeBounds );
			printf( "%d-th range: <%d %d %d> <%d %d %d>\n", rangeId,
									rangeBounds[0], rangeBounds[1], rangeBounds[2],
									rangeBounds[3], rangeBounds[4], rangeBounds[5] );

			rangeId++;
		}
	}// end function

};
#endif
/* FEL_ENCODER_H_ */
