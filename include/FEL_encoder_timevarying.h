/*
 * FEL_encoder_timevarying.h
 *
 *  Created on: Mar 18, 2012
 *      Author: abon
 */

#ifndef FEL_ENCODER_TIMEVARYING_H_
#define FEL_ENCODER_TIMEVARYING_H_

#include <list>
#include <pthread.h>

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
#include "FEL_spacetree_timevarying.h"
#include "FEL_segmenttree.h"
#include "FEL_matlab.h"

#define TMAX 100

template <class T>
class FEL_encoder_timevarying
{
	enum dataTypes { SCALARDATA = 0,
					 VECTORDATA = 1,
					 TIMEVARYINGSCALARDATA = 2,
					 TIMEVARYINGVECTORDATA = 3 };

	enum encodingTypes { REGULAR = 0,
						 QUADTREE = 1,
						 OCTREE = 2 };

	enum transformTypes { ROTATION = 0,
						  REFLECTION = 1 };

	enum domainCreationTypes{ SINGLERES = 0,
							  MULTIRES = 1
							};

	int nDim;
	float dataLow[3], dataHigh[3];
	int dataSize[3];

	int nTimeStep;

	int nTotalDomain;
	int nTotalRange, nLeafRange;

	float minDomainEntropy, maxDomainEntropy;
	bool isSharingNeighbor;

	ITL_histogram *histogram;

	list< FEL_range<T> > rangeBlockList;
	list< FEL_domain<T> > domainList;

	int globalRangeIndex;

	FEL_spacetree_timevarying<T> rangeBlockTree;

	ITL_partition_quadtree<T> quadtreePartitioner;
	ITL_partition_octree<T> octreePartitioner;

	FEL_segmenttree<T> domainTree;

	FEL_transformation trManager;

	Engine *engine;
	double errorList[1000];
	double errorColorList[3000];
	double matchFreqList[1000];
	double matchedDomainIdList[1000];
	double rangeLimitList[7000];
	double domainLimitList[7000];
	double domainEntropyList[1000];
	double queryDist[720];
	double resultDist[18000];			// 25 * 720
	double queryRegionLimitList[7];

public:

	int dataType;

	int smallestDomainSize[4];
	int smallestRangeSize[4];


	list<ITL_distribution> salientRangeDistributionList;

	bool isDomainPruningOn;

	//double resultRegionLimitList[175]; 	// 25 * 7
	float resultRegionLimitList[175]; 	// 25 * 7

	// Histogram properties
	int histMappingType;
	int nBin;

	// Domain creation parameters
	int domainCreationType;

	// Encoding parameters
	int encodingType;
	int distanceMeasureType;
	float matchErrorThreshold;
	float entropyThreshold;
	int nSegment;

	double domainCreationTime, domainOrgTime, encodeTime, searchTime;
	clock_t starttime, endtime;

	// Plotting
	bool matlab_plotting_on;
	FEL_matlab* plotter;

	ITL_field_regular<T>** inputDataField;
	ITL_field_regular<int>* inputBinField;

	//ITL_field_regular<T>* downsampledDataField;
	//ITL_field_regular<int>* downsampledBinField;
	//ITL_field_regular<T>** downsampledDataFieldArray;
	//ITL_field_regular<int>** downsampledBinFieldArray;

	// Constructor
	FEL_encoder_timevarying()
	{
		// Set default values to some parameters
		dataType = TIMEVARYINGSCALARDATA;
		nTimeStep = 1;
		encodingType = REGULAR;
		histMappingType = 2;
		domainCreationType = MULTIRES;
		distanceMeasureType = 0;
		isDomainPruningOn = false;

		// Initialize pointers to NULL
		inputDataField = NULL;
		inputBinField = NULL;
		histogram = NULL;

	}// end constructor

	FEL_encoder_timevarying( int datatype, int ndim, int *domainsize, int *rangesize, int nt )
	{
		// Get data type
		dataType = datatype;

		// Get data properties
		nTimeStep = nt;
		nDim = ndim;


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

		// Initialize pointers to NULL
		inputDataField = NULL;
		inputBinField = NULL;
		histogram = NULL;


	}// End constructor

	void
	initialize( int datatype, int ndim, int *domainsize, int *rangesize, int nt )
	{
		// Get data type
		dataType = datatype;

		// Get data
		nTimeStep = nt;
		nDim = ndim;

		// Record smallest domain and range sizes
		for( int i=0; i<nDim; i++ )
		{
			smallestDomainSize[i] = domainsize[i];
			smallestRangeSize[i] = rangesize[i];
		}

		// Initialize transformation class
		//trManager.initTransformationMatrices( 16 );

		// Initialize some other parameters
		nTotalDomain = 0;
		globalRangeIndex = 0;
		isSharingNeighbor = false;
		memset( matchFreqList, 0, 1000 );

	}// End constructor

	void
	setHistogram( ITL_histogram *hist, int nbin, int mappingtype )
	{
		histogram = hist;
		nBin = nbin;
		histMappingType = mappingtype;
	}

	void createDomains( ITL_field_regular<T>* inputDataField, int t = 0 )
	{

		fprintf( stderr, "Creating domains from time step %d ... \n", t );

		// Create/Update domains
		// Only multi-size algo is implemented for the time being
		createDomainsMultiSize( &inputDataField, t );


		// Update total number of domains
		nTotalDomain = domainList.size();
		fprintf( stderr, "%s: %d: %d domains created ...\n", __FILE__, __LINE__, nTotalDomain );
		cout << domainList.size() << endl;

	}

	void createDomainsMultiSize( ITL_field_regular<T>** inputDataField, int t = 0 )
	{
		float lowF[4];
		float highF[4];
		int dim[4], dsfDim[4];
		int dsfLowInt[4];
		int dsfHighInt[4];
		int resampledFieldSize[4];
		T* resampledData = NULL;
		ITL_field_regular<T>* downsampledDataField = new ITL_field_regular<T>();
		ITL_field_regular<int>* downsampledBinField = new ITL_field_regular<int>();

		fprintf( stderr, "Creating domains of different sizes from single resolution ... \n");

		// Down sample data to one level
		fprintf( stderr, "Down sampling field ... \n");
		(*inputDataField)->getSize( dataSize );
		FEL_util<T>::resampleField( (*inputDataField)->getDataFull(), dataType, dataSize,
									&resampledData, resampledFieldSize, 0.5f );

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
		cout << nDim << " " << nBin << " " << dsfDim[0] << " " << dsfDim[1] << " " << dsfDim[2] << endl;
		assert( histogram != NULL );
		if( histMappingType == 2 )
			FEL_util<T>::computeHistogramBinField2D( dataType,
													 &downsampledDataField, &downsampledBinField,
													 histogram, nBin,
													 dsfDim );
		else
			FEL_util<T>::computeHistogramBinField( dataType,
												   &downsampledDataField, &downsampledBinField,
												   histogram, nBin,
												   dsfDim );

		int m = ITL_util<int>::Min( downsampledBinField->getDataFull(), downsampledBinField->getSize() );
		int M = ITL_util<int>::Max( downsampledBinField->getDataFull(), downsampledBinField->getSize() );
		cout << "xx: " << m << " " << M << endl;

		// Delete down sampled data
		// From now on use only bin field
		delete downsampledDataField;

		// Create domains if first time step,
		// Otherwise just update the domains
		if( t == 0 )
		{
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
				fprintf( stderr, "Creating domains of size <%d, %d, %d> ...\n",
						 	 	 domainSize[0], domainSize[1], domainSize[2] );
				downsampledBinField->getBounds( dsfLowInt, dsfHighInt );
				createRegularDomains( dsfLowInt, dsfHighInt,
									  resampledFieldSize,
									  domainSize, i, t );

			}// end for : size i

		}// end if: t

		// Compute and update distributions of domains
		fprintf( stderr, "Updating domain distributions for time %d ... \n", t );
		updateDomainDistributions( &downsampledBinField, t );

		// Free temporary resources
		delete downsampledBinField;

	}// end function

	void createRegularDomains( int* dsfLowInt, int* dsfHighInt,
							   int* dsfDim,
							   int* domainSize,
							   int sizeIndex, int t )
	{
		int nDomain[nDim];
		float lowSub[nDim];
		float highSub[nDim];

		// Compute the number of domains
		nDomain[0] = (int)ceil( dsfDim[0]/(float)domainSize[0] );
		nDomain[1] = (int)ceil( dsfDim[1]/(float)domainSize[1] );
		nDomain[2] = (int)ceil( dsfDim[2]/(float)domainSize[2] );

		//#ifdef DEBUG_MODE
		printf( "Dimension of the downsampled field: %d, %d, %d\n", dsfDim[0], dsfDim[1], dsfDim[2] );
		printf( "Number of domains to be created: %d, %d, %d\n", nDomain[0], nDomain[1], nDomain[2] );
		printf( "Size of a domain: %d, %d, %d\n", domainSize[0], domainSize[1], domainSize[2] );
		//#endif

		// Nested loop for partitioning
		int domainIndex = 0;
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
					FEL_domain<T> newDomain( nDim, lowSub, highSub, sizeIndex, nBin, nTimeStep );

					// Store domain in list
					domainList.push_back( newDomain );

					// Increment the index to the next domain
					domainIndex ++;

				}
			}
		}

	}

	void
	updateDomainDistributions( ITL_field_regular<int>** downsampledBinField, int t = 0 )
	{
		int domainBounds[6];
		float lowSub[4];
		float highSub[4];
		int nPointDomain = 1;
		double* domainFreqList = new double[nBin];
		//double* domainTvFreqList = new double[nBin*nTimeStep];
		int* domainBinData = NULL;
		float domainEntropy;
		int nb, nt;

		int domainID = 0;
		for( typename list< FEL_domain<T> >::iterator domainIter = domainList.begin();
			 domainIter != domainList.end();
			 domainIter++ )
		{
			// Get domain block bounds and size
			domainIter->getBounds( lowSub, highSub );
			nPointDomain = domainIter->getSize();
			///cout << "1: " << domainID << " " << nPointDomain << endl;
			//cout << lowSub[0] << " " << lowSub[1] << " " << lowSub[2] << endl;
			//cout << highSub[0] << " " << highSub[1] << " " << highSub[2] << endl;

			// Get domain bin data
			domainBinData = new int[nPointDomain];
			(*downsampledBinField)->getDataBetween( lowSub, highSub, domainBinData );
			int m = ITL_util<int>::Min( domainBinData, nPointDomain );
			int M = ITL_util<int>::Max( domainBinData, nPointDomain );
			//cout << "2: " << m << " " << M << endl;

			// Convert to distribution
			FEL_util<T>::computeHistogramFromBinField( domainBinData, nPointDomain, nBin, domainFreqList );
			//cout << "3: " << nPointDomain << endl;

			// Get domain distribution
			//domainIter->getDistribution( domainTvFreqList, &nb, &nt );
			//cout << "4: " << nPointDomain << endl;

			// Update distribution for particular time step
			//memcpy( domainTvFreqList + t*nBin, domainFreqList, sizeof(double)*nBin );
			//cout << "5: " << nPointDomain << endl;

			// Store updated domain distribution
			assert( domainFreqList != NULL );
			domainIter->setDistribution( domainFreqList, nBin, t );
			//cout << "6: " << nPointDomain << endl;

			domainID ++;

			// Free temporary resources
			delete [] domainBinData;

		}// end for : list of domains

		delete [] domainFreqList;
		//delete [] domainTvFreqList;

	}// end function

	void computeDomainEntropies()
	{
		int domainBounds[6];
		int nPointDomain = 1;
		float timeEntropyList[nTimeStep];
		double domainTvFreqList[nBin*nTimeStep];
		float domainEntropy;
		int nb, nt;

		int domainID = 0;
		for( typename list< FEL_domain<T> >::iterator domainIter = domainList.begin();
			 domainIter != domainList.end();
			 domainIter++ )
		{
			// Get domain block bounds and distribution
			//cout << "1" << endl;
			domainIter->getBounds( domainBounds );
			cout << "2" << endl;
			domainIter->getTimeVaryingDistribution( domainTvFreqList, &nb, &nt );

			// Compute and store domain entropy
			cout << "3" << endl;
			domainEntropy = computeEntropy_TimeHistogram( domainTvFreqList, nBin, nTimeStep, false, timeEntropyList );
			cout << "3.5" << endl;
			domainIter->setEntropy( domainEntropy );
			cout << "4" << endl;
			domainIter->setTimevaryingEntropy( timeEntropyList );

			cout << "4.5" << endl;
			if( domainID == 0 )
			{
				minDomainEntropy = maxDomainEntropy = domainEntropy;
			}
			else
			{
				if( domainEntropy < minDomainEntropy ) minDomainEntropy = domainEntropy;
				if( domainEntropy > maxDomainEntropy ) maxDomainEntropy = domainEntropy;
			}

			// Store domain entropies for display
			//domainEntropyList[domainID] = domainEntropy;
			domainID ++;

			// Free temporary resources
			//delete [] domainBinData;
		}

	}// end function

	float
	computeEntropy_TimeHistogram( double* tvFreqList, int nBin, int ntime, bool isNormalize, float* timeEntropyList )
	{
		float sum = 0;
		for( int t=0; t<ntime; t++ )
		{
			timeEntropyList[t] = (float)ITL_entropycore::computeEntropy_HistogramBased2( tvFreqList + t*nBin, nBin, isNormalize );
			sum += timeEntropyList[t];
		}

		return ( sum / (float)ntime );

	}// end function

	float
	computeTimeWiseEntropyDifferenceMax( float* hList1, float* hList2, int ntime )
	{
		float max = 0;
		float diff;
		for( int t=0; t<ntime; t++ )
		{
			diff = abs( hList1[t] - hList2[t] );
			if( diff > max )
				max = diff;
		}

		return max;

	}// end function


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
		//float entropyThreshold = 0.01f * domainEntropyRange;
		float entropyThreshold = 0.05f * domainEntropyRange;

		//float lastEntropy, curEntropy;
		float lastEntropyList[nTimeStep], curEntropyList[nTimeStep];

		for( typename list< FEL_domain<T> >::iterator domainIter = domainList.begin();
			 domainIter != domainList.end();
			 domainIter++ )
		{
			// Get entropy list of current domain
			domainIter->getTimevaryingEntropy( curEntropyList );

			if( i==0 )
			{
				//lastEntropy = curEntropy;
				memcpy( lastEntropyList, curEntropyList, sizeof(float)*nTimeStep );
				reducedDomainList.push_back( *domainIter );
			}
			else
			{
				if( computeTimeWiseEntropyDifferenceMax( curEntropyList, lastEntropyList, nTimeStep ) > entropyThreshold )
				{
					//lastEntropy = curEntropy;
					memcpy( lastEntropyList, curEntropyList, sizeof(float)*nTimeStep );
					reducedDomainList.push_back( *domainIter );
				}
			}

			i++;

		}// end for

		// Update list of domains
		domainList.clear();
		domainList = reducedDomainList;

	}

	void organizeDomains()
	{
		starttime = ITL_util<float>::startTimer();

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
			pruneDomains_Entropybased();
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
	encode_computeDistribution( ITL_field_regular<T>** inputdatafield, int t = 0 )
	{
		// Shallow copy
		inputDataField = inputdatafield;

		starttime = ITL_util<float>::startTimer();
		if( encodingType == REGULAR )
		{
			fprintf( stderr, "Initializing encoding with regular partitioning ... \n" );
			encode_Regular_computeDistribution( t );
		}
		else
		{
			fprintf( stderr, "Initializing encoding with quadtree partitioning ... \n" );
			encode_Spacetree_computeDistribution( t );
		}
		encodeTime = ITL_util<float>::endTimer( starttime );
		fprintf( stderr, "Total encoding time: %g\n", encodeTime );

	}

	void
	encode_Regular_computeDistribution( int t = 0 )
	{
		int nRange[nDim];
		float lowF[3], highF[3];
		float lowSub[nDim];
		float highSub[nDim];
		int lowInt[nDim], highInt[nDim];
		int dim[nDim];
		(*inputDataField)->getBounds( lowInt, highInt );
		(*inputDataField)->getBounds( lowF, highF );
		(*inputDataField)->getSize( dim );

		//cout << dim[0] << " " << dim[1] << " " << dim[2] << endl;
		//cout << lowF[0] << " " << lowF[1] << " " << lowF[2] << endl;
		//cout << highF[0] << " " << highF[1] << " " << highF[2] << endl;

		// Convert to distribution
		inputBinField = new ITL_field_regular<int>( nDim, lowF, highF );
		if( histMappingType == 2 )
			FEL_util<T>::computeHistogramBinField2D( dataType,
													 &inputDataField, &inputBinField,
													 histogram, nBin, dim );
		else
			FEL_util<T>::computeHistogramBinField( dataType,
												   &inputDataField, &inputBinField,
												   histogram, nBin, dim );

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
		double rangeTvFreqList[nBin*nTimeStep];

		// Create the range blocks by subdividing the space
		// for the first time step
		// For the following time steps, only update distributions

		cout << "creating ranges / updating: " << t << endl;
		if( t == 0 )
		{
			for( int z=0; z<nRange[2]; z++ )
			{
				// Determine [zmin, zmax] for the next range
				if( dim[2] == 2 )
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
						FEL_util<int>::computeHistogramFromBinField( rangeBinData, nPointRange, nBin, rangeFreqList );
						delete [] rangeBinData;

						// Create a range block and
						// set properties of the range block including distribution
						FEL_range<T> nextRange( nDim, lowSub, highSub, nBin, nTimeStep );
						nextRange.setDistribution( rangeFreqList, nBin, t );

						//int nb, nt;
						//nextRange.getTimevaryingDistribution( rangeTvFreqList, &nb, &nt );
						//cout << "in create: " << nb << " " << nt << endl;

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

		}// end if : for t
		else
		{
			updateRangeDistributions( &inputBinField, t );
		}


		delete inputBinField;

	}// end function

	void
	updateRangeDistributions( ITL_field_regular<int>** inputBinField, int t = 0 )
	{
		int rangeBounds[6];
		float lowSub[4];
		float highSub[4];
		int nPointRange = 1;
		double rangeFreqList[nBin];
		double rangeTvFreqList[nBin*nTimeStep];
		int* rangeBinData = NULL;
		int nb, nt;

		int rangeID = 0;
		for( typename list< FEL_range<T> >::iterator rangeIter = rangeBlockList.begin();
			 rangeIter != rangeBlockList.end();
			 rangeIter++ )
		{
			// Get domain block bounds and size
			rangeIter->getBounds( lowSub, highSub );
			nPointRange = rangeIter->getSize();
			//cout << lowSub[0] << " " << lowSub[1] << " " << lowSub[2] << endl;
			//cout << highSub[0] << " " << highSub[1] << " " << highSub[2] << endl;
			//rangeIter->getTimevaryingDistribution( rangeTvFreqList, &nb, &nt );
			//cout << "in update: " << nb << " " << nt << endl;

			// Get domain bin data
			rangeBinData = new int[nPointRange];
			(*inputBinField)->getDataBetween( lowSub, highSub, rangeBinData );

			// Convert to distribution
			FEL_util<double>::computeHistogramFromBinField( rangeBinData, nPointRange, nBin, rangeFreqList );
			#ifdef DEBUG_MODE
			double m = ITL_util<double>::Min( rangeFreqList, nBin );
			double M = ITL_util<double>::Max( rangeFreqList, nBin );
			fprintf( stderr, "%d-th Range distribution min/max: %g, %g\n", rangeID, m, M );
			#endif


			// Get domain distribution
			//rangeIter->getTimeVaryingDistribution( rangeTvFreqList, &nb, &nt );

			// Update distribution for particular time step
			//memcpy( rangeTvFreqList + t*nBin, rangeFreqList, sizeof(double)*nBin );

			// Store updated domain distribution
			rangeIter->setDistribution( rangeFreqList, nBin, t );

			rangeID ++;

			// Free temporary resources
			delete [] rangeBinData;

		}// end for : list of domains

	}// end function

	void
	encode_Spacetree_computeDistribution( int t = 0  )
	{
		float lowF[3], highF[3];
		int dim[3];

		// Convert to distribution
		inputBinField = new ITL_field_regular<int>( nDim, lowF, highF );
		(*inputDataField)->getSize( dim );
		if( histMappingType == 2 )
			FEL_util<T>::computeHistogramBinField2D( dataType,
													 &inputDataField, &inputBinField,
													 histogram, nBin, dim );
		else
			FEL_util<T>::computeHistogramBinField( dataType,
												   &inputDataField, &inputBinField,
												   histogram, nBin, dim );

		// Initialize range block tree
		rangeBlockTree.initialize( &inputBinField,
								   &domainList,
							 	   &histogram, nBin,
								   smallestRangeSize,
								   distanceMeasureType,
								   encodingType,
								   matchErrorThreshold,
								   entropyThreshold,
								   &trManager );

		// Create root node (represents entire field)
		fprintf( stderr, "Creating range blobk quadtree root ... \n");
		inputBinField->getBounds( lowF, highF );
		rangeBlockTree.initTree( lowF, highF );

		// Recursively create range blocks
		fprintf( stderr, "Recursively creating range block quadtree ... \n");
		rangeBlockTree.createTree( rangeBlockTree.getRoot(), lowF, highF, 0 );

		// Recursively compute distributions of the nodes
		rangeBlockTree.fillTree( rangeBlockTree.getRoot(), 0 );

		// Obtain total number of ranges
		nTotalRange = rangeBlockTree.getNumRangeBlock();
		nLeafRange = rangeBlockTree.getNumLeafRangeBlock();
		fprintf( stderr, "%s: %d: %d rangeblocks created ...\n", __FILE__, __LINE__, nTotalRange );
		fprintf( stderr, "%s: %d: %d leaf nodes created ...\n", __FILE__, __LINE__, nLeafRange );

		delete inputBinField;

	}

	void
	encode_matchDistribution( int t = 0 )
	{
		float lowF[nDim], highF[nDim];
		int dim[nDim];

		if( encodingType == REGULAR )
		{
			//fprintf( stderr, "Initializing encoding with regular partitioning ... \n" );
			//encode_Regular_MatchDistribution();

			fprintf( stderr, "Initializing encoding with regular partitioning and partial matching ... \n" );
			encode_Regular_MatchDistribution_Part();
		}
		else
		{
			fprintf( stderr, "Initializing encoding with quadtree partitioning ... \n" );
			encode_Spacetree_MatchDistribution();
		}

	}

	void
	encode_Regular_MatchDistribution()
	{
		// Iterate though each range
		// and find the nearest domain
		int nb, nt;
		int bestMatch = -1;
		double matchError = -1;
		int rangeID = 0;
		int optRot = 0;
		bool optRef = false;
		float rangeLow[3];
		float rangeHigh[3];
		double* rangeTvFreqList = new double[nBin*nTimeStep];

		for( typename list< FEL_range<T> >::iterator rangeBlockIter = rangeBlockList.begin();
			 rangeBlockIter != rangeBlockList.end();
			 rangeBlockIter++ )
		{
			fprintf( stderr, "Finding best match for the next range ... \n" );

			//rangeBlockIter->getBounds( rangeLow, rangeHigh );
			//cout << rangeLow[0] << " " << rangeLow[1] << " " << rangeLow[2] << endl;
			//cout << rangeHigh[0] << " " << rangeHigh[1] << " " << rangeHigh[2] << endl;

			// Get Range distribution
			rangeBlockIter->getTimevaryingDistribution( rangeTvFreqList, &nb, &nt );
			//cout << nb << " " << nt << endl;
			#ifdef DEBUG_MODE
			double m = ITL_util<double>::Min( rangeTvFreqList, nBin*nTimeStep );
			double M = ITL_util<double>::Max( rangeTvFreqList, nBin*nTimeStep );
			fprintf( stderr, "%s: %d: %d-th Range distribution min/max: %g, %g\n", __FILE__, __LINE__, rangeID, m, M );
			#endif

			//bestMatch = matchRangeToDomains( rangeBlockIter, &matchError, &optimalTrMatrix );

			optRot = 0;
			optRef = false;
			bestMatch = matchRangeToDomainsByTimeHistogram( rangeTvFreqList, &matchError, &optRot, &optRef );

			//#ifdef DEBUG_MODE
			fprintf( stderr, "%d-th range maps to %d-th domain with %g <%d, %d> error ...\n",
					 rangeID, bestMatch, matchError, optRot, (int)optRef );
			//#endif

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
				plotter->errorList[rangeID] = matchError;
				plotter->matchedDomainIdList[rangeID] = bestMatch;
				plotter->domainMatchFreqList[bestMatch] ++;
			}

			// Increment range ID
			rangeID ++;
		}


		for( typename list< FEL_domain<T> >::iterator domIter = domainList.begin();
			 domIter != domainList.end();
			 domIter++ )
		{
			domIter->sortCodeBookEntries();
		}

		delete [] rangeTvFreqList;

	}

	void
	encode_Spacetree_MatchDistribution()
	{
		float lowF[3], highF[3];

		// Recursively compute distributions of the nodes
		rangeBlockTree.matchTree( rangeBlockTree.getRoot() );

	}

	int
	matchRangeToDomainsByTimeHistogram( double* rangeTvFreqList, double *error, int* rotateAmount, bool* isReflected )
	{
		int nb, nt;
		int bestMatchedDomainId = -1;
		double e = 100;
		double minError = 100;
		double domainTvFreqList[nBin*nTimeStep];
		int domainID = 0;
		int rotAmount = 0;
		bool isRef = false;

		cout << "Number of domains " << domainList.size() << endl;
		for( typename list< FEL_domain<T> >::iterator domainIter = domainList.begin();
			 domainIter != domainList.end();
			 domainIter++ )
		{
			domainIter->getTimeVaryingDistribution( domainTvFreqList, &nb, &nt );
			#ifdef DEBUG_MODE
			double m = ITL_util<double>::Min( domainTvFreqList, nBin*nTimeStep );
			double M = ITL_util<double>::Max( domainTvFreqList, nBin*nTimeStep );
			fprintf( stderr, "%d-th Domain distribution min/max: %g, %g\n", domainID, m, M );
			#endif

			e =  matchRangeToDomainByTimeHistogram( rangeTvFreqList, domainTvFreqList, &rotAmount, &isRef );
			if( e < minError )
			{
				minError = e;
				(*rotateAmount) = rotAmount;
				(*isReflected) = isRef;
				bestMatchedDomainId = domainID;
			}

			domainID ++;
		}

		*error = minError;
		return bestMatchedDomainId;

	}

	double
	matchRangeToDomainByTimeHistogram( double* rangeTvFreqList, double* domainTvFreqList,
									   int* rotateAmount, bool* isReflected )
	{
		double* domainTvFreqListCopy = new double[nBin*nTimeStep];
		double* domainTvFreqListRev = new double[nBin*nTimeStep];
		double* buffer = new double[nTimeStep];
		double* buffer2 = new double[(nBin-1)*nTimeStep];

		// Apply various transformations to the domain
		// Compute distributions of the transformed domain
		// and compare to the range distribution
		double matchErrorPrcnt = 100;
		double minErrorPrcnt = 100;
		int optRot = 0;
		bool optRef = false;

		memcpy( domainTvFreqListCopy, domainTvFreqList, sizeof(double)*(nBin*nTimeStep) );
		for( int iT = 0; iT<nBin; iT++ )
		{
			// Compare with range
			matchErrorPrcnt = FEL_util<double>::computeTimeHistogramMatchError( distanceMeasureType,
																				rangeTvFreqList, domainTvFreqListCopy,
																				nBin, nTimeStep );

			if( matchErrorPrcnt < minErrorPrcnt )
			{
				minErrorPrcnt = matchErrorPrcnt;
				//*optTransformationMatrix = trManager.rotMatrixArray[iT];
				optRot = iT;
				optRef = false;
			}

			// Place right(top)most to buffer
			for( int t=0; t<nTimeStep; t++ )
				buffer[t] = domainTvFreqListCopy[t*nBin+(nBin-1)];

			// Move everybody else by 1 bit
			//memcpy( buffer2, domainFreqListCopy, sizeof(float)*(nBin-1) );
			//memcpy( domainFreqListCopy+1, buffer2, sizeof(float)*(nBin-1) );
			for( int t=0; t<nTimeStep; t++ )
			{
				memcpy( buffer2+t*(nBin-1), domainTvFreqListCopy+t*nBin, sizeof(double)*(nBin-1) );
				memcpy( domainTvFreqListCopy+t*nBin+1, buffer2+t*(nBin-1), sizeof(double)*(nBin-1) );
			}
			//for( int j=nBin-1; j>=1; j-- )
			//	domainFreqListCopy[j] = domainFreqListCopy[j-1];

			// Fill back the leftmost
			for( int t=0; t<nTimeStep; t++ )
				domainTvFreqListCopy[t*nBin] = buffer[t];

			//printf( "\n" );
			//for( int k=0; k<nBin; k++ )
			//	printf( "%g, ", domainFreqListCopy[k] );
			//printf( "\n" );
		}

		// Reflect distribution
		for( int t=0; t<nTimeStep; t++ )
		{
			for( int i=0; i<nBin; i++ )
				domainTvFreqListRev[t*nBin+i] = domainTvFreqList[t*nBin+(nBin-i-1)];
		}
		memcpy( domainTvFreqListCopy, domainTvFreqListRev, sizeof(double)*(nBin*nTimeStep) );

		//printf( "Ref: \n" );
		//for( int k=0; k<nBin; k++ )
		//	printf( "%g, ", domainFreqListCopy[k] );
		//printf( "\n" );

		for( int iT = 0; iT<nBin; iT++ )
		{
			// Compare with range distribution
			matchErrorPrcnt = FEL_util<double>::computeTimeHistogramMatchError( distanceMeasureType,
																				rangeTvFreqList, domainTvFreqListCopy,
																				nBin, nTimeStep );

			if( matchErrorPrcnt < minErrorPrcnt )
			{
				minErrorPrcnt = matchErrorPrcnt;
				optRot = iT;
				optRef = true;
			}

			// Place right(top)most to buffer
			for( int t=0; t<nTimeStep; t++ )
				buffer[t] = domainTvFreqListCopy[t*nBin+(nBin-1)];

			// Move everybody else by 1 bit
			//memcpy( buffer2, domainFreqListCopy, sizeof(float)*(nBin-1) );
			//memcpy( domainFreqListCopy+1, buffer2, sizeof(float)*(nBin-1) );
			for( int t=0; t<nTimeStep; t++ )
			{
				memcpy( buffer2+t*(nBin-1), domainTvFreqListCopy+t*nBin, sizeof(double)*(nBin-1) );
				memcpy( domainTvFreqListCopy+t*nBin+1, buffer2+t*(nBin-1), sizeof(double)*(nBin-1) );
			}
			//for( int j=nBin-1; j>=1; j-- )
			//	domainFreqListCopy[j] = domainFreqListCopy[j-1];

			// Fill back the leftmost
			for( int t=0; t<nTimeStep; t++ )
				domainTvFreqListCopy[t*nBin] = buffer[t];

			//printf( "\n" );
			//for( int k=0; k<nBin; k++ )
			//	printf( "%g, ", domainFreqListCopy[k] );
			//printf( "\n" );
		}

		(*rotateAmount) = optRot;
		(*isReflected) = optRef;

		// Free temporary resources
		delete [] domainTvFreqListCopy;
		delete [] domainTvFreqListRev;
		delete [] buffer;
		delete [] buffer2;

		return minErrorPrcnt;

	}// End function

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	encode_Regular_MatchDistribution_Part()
	{
		// Iterate though each range
		// and find the nearest domain
		int nb, nt;
		int bestMatch = -1;
		double matchError = -1;
		int rangeID = 0;
		int optRot = 0;
		bool optRef = false;
		float rangeLow[3];
		float rangeHigh[3];
		double* rangeTvFreqList = new double[nBin*nTimeStep];
		int bestMatchedDomainIDList[nSegment];
		int bestMatchedSegmentIDList[nSegment];
		double matchErrorList[nSegment];
		int rotAmountList[nSegment];
		bool refFlagList[nSegment];

		for( typename list< FEL_range<T> >::iterator rangeBlockIter = rangeBlockList.begin();
			 rangeBlockIter != rangeBlockList.end();
			 rangeBlockIter++ )
		{
			fprintf( stderr, "Finding best match for the next range ... \n" );

			// Get Range distribution
			rangeBlockIter->getTimevaryingDistribution( rangeTvFreqList, &nb, &nt );
			#ifdef DEBUG_MODE
			double m = ITL_util<double>::Min( rangeTvFreqList, nBin*nTimeStep );
			double M = ITL_util<double>::Max( rangeTvFreqList, nBin*nTimeStep );
			fprintf( stderr, "%s: %d: %d-th Range distribution min/max: %g, %g\n", __FILE__, __LINE__, rangeID, m, M );
			#endif

			// Match range with next domain
			matchRangeToDomainsByTimeHistogramPart( nSegment,
													rangeTvFreqList,
													matchErrorList,
													bestMatchedDomainIDList,
													bestMatchedSegmentIDList,
													rotAmountList,
													refFlagList );
			//#ifdef DEBUG_MODE
			fprintf( stderr, "Matching information for %d-th range:\n", rangeID );
			for( int iS=0; iS<nSegment; iS++ )
				fprintf( stderr, "<%d, %d><%d, %d>%g\t", bestMatchedDomainIDList[iS],
														 bestMatchedSegmentIDList[iS],
														 rotAmountList[iS],
														 (int)refFlagList[iS],
														 matchErrorList[iS] );
			fprintf( stderr, "\n" );
			//#endif

			// Add range information to domain
			fprintf( stderr, "Adding range information to the matched domain segments (possibly from different domains) ... \n" );
			for( int iS=0; iS<nSegment; iS++ )
			{
				typename list<FEL_domain<T> >::iterator tmpIter = domainList.begin();
				for( int iD=0; iD<bestMatchedDomainIDList[iS]; iD++ )
					tmpIter++;
				tmpIter->addMatchingRangeSegment( rangeID, iS,
												  bestMatchedSegmentIDList[iS],
												  matchErrorList[iS],
												  rotAmountList[iS],
												  refFlagList[iS] );
			}

			// Add domain information to range
			fprintf( stderr, "Adding domain informations to the matched range ... \n" );
			rangeBlockIter->setNumSegment( nSegment );
			cout << "22" << endl;
			rangeBlockIter->setMatchedDomainSegmentId( bestMatchedDomainIDList, bestMatchedSegmentIDList );
			cout << "23" << endl;
			rangeBlockIter->setMatchError( matchErrorList  );
			cout << "24" << endl;
			rangeBlockIter->setMatchTransformation( rotAmountList, refFlagList );
			cout << "25" << endl;

			// Store Matching information
			// for Matlab plotting
			errorList[rangeID] = matchError;
			errorColorList[rangeID*3+2] = 0.3 + 0.7*(matchError / 100);
			errorColorList[rangeID*3] = errorColorList[rangeID*3+1] = 0.005;
			matchedDomainIdList[rangeID] = bestMatch;
			matchFreqList[bestMatch] ++;

			// Increment range ID
			rangeID ++;
		}


		for( typename list< FEL_domain<T> >::iterator domIter = domainList.begin();
			 domIter != domainList.end();
			 domIter++ )
		{
			domIter->sortCodeBookEntries();
		}

		delete [] rangeTvFreqList;

	}



	void
	matchRangeToDomainsByTimeHistogramPart( int nSegment,
											double* rangeTvFreqList,
											double *error,
											int* matchingDomain,
											int* matchingSegment,
											int* rotateAmount,
											bool* isReflected )
	{
		int nb, nt;
		double domainTvFreqList[nBin*nTimeStep];

		int bestMatchedDomainId[nSegment];
		int bestMatchedSegmentId[nSegment];
		//double matchErrorList[nSegment];
		double minErrorList[nSegment];
		int rotAmount[nSegment];
		bool refList[nSegment];
		int domainID;
		memset( bestMatchedDomainId, -1, sizeof(int)*nSegment );
		//memset( matchErrorList, 100, sizeof(double)*nSegment );
		memset( minErrorList, 100, sizeof(double)*nSegment );

		int segmentIdList[nSegment];
		double errorList[nSegment];
		int rotAmountList[nSegment];
		bool refFlagList[nSegment];
		memset( errorList, 100, sizeof(double)*nSegment );
		cout << "Number of domains " << domainList.size() << endl;

//#ifdef USING_PTHREAD
		int nIter = domainList.size();

//#else
		domainID = 0;
		for( typename list< FEL_domain<T> >::iterator domainIter = domainList.begin();
			 domainIter != domainList.end();
			 domainIter++ )
		{
			domainIter->getTimeVaryingDistribution( domainTvFreqList, &nb, &nt );
			//#ifdef DEBUG_MODE
			double m = ITL_util<double>::Min( domainTvFreqList, nBin*nTimeStep );
			double M = ITL_util<double>::Max( domainTvFreqList, nBin*nTimeStep );
			fprintf( stderr, "%d-th Domain distribution min/max: %g, %g\n", domainID, m, M );
			//#endif

			matchRangeToDomainByTimeHistogramPart( nSegment,
												   rangeTvFreqList, domainTvFreqList,
												   errorList, segmentIdList,
												   rotAmountList, refFlagList );

			// For each range segment, update the best match information
			// if found
			for( int p=0; p<nSegment; p++ )
			{
				if( errorList[p] < minErrorList[p] )
				{
					minErrorList[p] = errorList[p];
					bestMatchedDomainId[p] = domainID;
					bestMatchedSegmentId[p] = segmentIdList[p];
					rotAmount[p] = rotAmountList[p];
					refList[p] = refFlagList[p];
				}
			}

			domainID ++;

		}// end for : domain list iterator

//#endif

		memcpy( error, minErrorList, sizeof(double)*nSegment );
		memcpy( rotateAmount, rotAmount, sizeof(int)*nSegment );
		memcpy( isReflected, refList, sizeof(bool)*nSegment );
		memcpy( matchingDomain, bestMatchedDomainId, sizeof(int)*nSegment );
		memcpy( matchingSegment, bestMatchedSegmentId, sizeof(int)*nSegment );

	}

	void
	matchRangeToDomainByTimeHistogramPart( int nSegment,
										   double* rangeTvFreqList, double* domainTvFreqList,
										   double* matchError, int* matchingSegmentId,
									   	   int* rotateAmount, bool* isReflected )
	{
		// Number of time steps in each time segment
		int ntPerSegment = nTimeStep / nSegment;

		double errorPerSegment[nSegment];
		int rotateAmountPerSegment[nSegment];
		int isReflectedPerSegment[nSegment];
		int bestmatchingSegmentID[nSegment];

		// Find best match for each segment of current range block
		// Each range segment is compared against each domain segment
		for( int p=0; p<nSegment; p++ )
		{
			int tStartSegmentRange = p*ntPerSegment;

			// Find which domain segment is best for this
			// range segment
			float minError2 = 100;
			int bestMatchSegmentId;
			int rot;
			bool ref;
			float error;
			for( int q=0; q<nSegment; q++ )
			{
				int tStartSegmentDomain = q*ntPerSegment;
				float error = matchRangeToDomainByTimeHistogram( rangeTvFreqList+nBin*tStartSegmentRange,
															     domainTvFreqList+nBin*tStartSegmentDomain,
															     &rot, &ref,
															     ntPerSegment );

				// Update record for current range segment
				if( error < minError2 )
				{
					minError2 = error;
					errorPerSegment[p] = error;
					rotateAmountPerSegment[p] = rot;
					isReflectedPerSegment[p] = ref;
					bestmatchingSegmentID[p] = q;
				}

			}// end inner for

		}// end outer for

		memcpy( matchError, errorPerSegment, sizeof(double)*nSegment );
		memcpy( rotateAmount, rotateAmountPerSegment, sizeof(int)*nSegment );
		memcpy( isReflected, isReflectedPerSegment, sizeof(bool)*nSegment );
		memcpy( matchingSegmentId, bestmatchingSegmentID, sizeof(int)*nSegment );

	}// End function

	/**
	 * Special: uses own ntime
	 */
	double
	matchRangeToDomainByTimeHistogram( double* rangeTvFreqList, double* domainTvFreqList,
									   int* rotateAmount, bool* isReflected, int ntime )
	{
		double* domainTvFreqListCopy = new double[nBin*ntime];
		double* domainTvFreqListRev = new double[nBin*ntime];
		double* buffer = new double[ntime];
		double* buffer2 = new double[(nBin-1)*ntime];

		// Apply various transformations to the domain
		// Compute distributions of the transformed domain
		// and compare to the range distribution
		double matchErrorPrcnt = 100;
		double minErrorPrcnt = 100;
		int optRot = 0;
		bool optRef = false;

		memcpy( domainTvFreqListCopy, domainTvFreqList, sizeof(double)*(nBin*ntime) );
		for( int iT = 0; iT<nBin; iT++ )
		{
			// Compare with range
			matchErrorPrcnt = FEL_util<double>::computeTimeHistogramMatchError( distanceMeasureType,
																				rangeTvFreqList, domainTvFreqListCopy,
																				nBin, ntime );

			if( matchErrorPrcnt < minErrorPrcnt )
			{
				minErrorPrcnt = matchErrorPrcnt;
				//*optTransformationMatrix = trManager.rotMatrixArray[iT];
				optRot = iT;
				optRef = false;
			}

			// Place right(top)most to buffer
			for( int t=0; t<ntime; t++ )
				buffer[t] = domainTvFreqListCopy[t*nBin+(nBin-1)];

			// Move everybody else by 1 bit
			//memcpy( buffer2, domainFreqListCopy, sizeof(float)*(nBin-1) );
			//memcpy( domainFreqListCopy+1, buffer2, sizeof(float)*(nBin-1) );
			for( int t=0; t<ntime; t++ )
			{
				memcpy( buffer2+t*(nBin-1), domainTvFreqListCopy+t*nBin, sizeof(double)*(nBin-1) );
				memcpy( domainTvFreqListCopy+t*nBin+1, buffer2+t*(nBin-1), sizeof(double)*(nBin-1) );
			}
			//for( int j=nBin-1; j>=1; j-- )
			//	domainFreqListCopy[j] = domainFreqListCopy[j-1];

			// Fill back the leftmost
			for( int t=0; t<ntime; t++ )
				domainTvFreqListCopy[t*nBin] = buffer[t];

			//printf( "\n" );
			//for( int k=0; k<nBin; k++ )
			//	printf( "%g, ", domainFreqListCopy[k] );
			//printf( "\n" );
		}

		// Reflect distribution
		for( int t=0; t<ntime; t++ )
		{
			for( int i=0; i<nBin; i++ )
				domainTvFreqListRev[t*nBin+i] = domainTvFreqList[t*nBin+(nBin-i-1)];
		}
		memcpy( domainTvFreqListCopy, domainTvFreqListRev, sizeof(double)*(nBin*ntime) );

		//printf( "Ref: \n" );
		//for( int k=0; k<nBin; k++ )
		//	printf( "%g, ", domainFreqListCopy[k] );
		//printf( "\n" );

		for( int iT = 0; iT<nBin; iT++ )
		{
			// Compare with range distribution
			matchErrorPrcnt = FEL_util<double>::computeTimeHistogramMatchError( distanceMeasureType,
																				rangeTvFreqList, domainTvFreqListCopy,
																				nBin, ntime );

			if( matchErrorPrcnt < minErrorPrcnt )
			{
				minErrorPrcnt = matchErrorPrcnt;
				optRot = iT;
				optRef = true;
			}

			// Place right(top)most to buffer
			for( int t=0; t<ntime; t++ )
				buffer[t] = domainTvFreqListCopy[t*nBin+(nBin-1)];

			// Move everybody else by 1 bit
			//memcpy( buffer2, domainFreqListCopy, sizeof(float)*(nBin-1) );
			//memcpy( domainFreqListCopy+1, buffer2, sizeof(float)*(nBin-1) );
			for( int t=0; t<ntime; t++ )
			{
				memcpy( buffer2+t*(nBin-1), domainTvFreqListCopy+t*nBin, sizeof(double)*(nBin-1) );
				memcpy( domainTvFreqListCopy+t*nBin+1, buffer2+t*(nBin-1), sizeof(double)*(nBin-1) );
			}
			//for( int j=nBin-1; j>=1; j-- )
			//	domainFreqListCopy[j] = domainFreqListCopy[j-1];

			// Fill back the leftmost
			for( int t=0; t<ntime; t++ )
				domainTvFreqListCopy[t*nBin] = buffer[t];

			//printf( "\n" );
			//for( int k=0; k<nBin; k++ )
			//	printf( "%g, ", domainFreqListCopy[k] );
			//printf( "\n" );
		}

		(*rotateAmount) = optRot;
		(*isReflected) = optRef;

		// Free temporary resources
		delete [] domainTvFreqListCopy;
		delete [] domainTvFreqListRev;
		delete [] buffer;
		delete [] buffer2;

		return minErrorPrcnt;

	}// End function


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
		FEL_range<T> returnRange;
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
		for( int t=0; t<nTimeStep; t++ )
		{
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
							nextRange = searchRangeByID( rangeIdList[iR], t );
							assert( nextRange != NULL );
							//cout << "here2" << endl;

						}
						else
						{
							nextRange = rangeBlockTree[t].searchTreeByRangeID( rangeBlockTree[t].getRoot(), rangeIdList[iR] );
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
							//returnList->push_back( returnRange );
							returnList.push_back( returnRange );

							nReturnedRange ++;

						}// end if

						//delete [] rangeIdList;

					}// end if-else

				}// end for : iR

			}// end for : domain iterator

		}// end for : t

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
		cout << "c" << endl;

		// For each returned domain, find the corresponding ranges
		nReturnedRange = 0;
		for( int t=0; t<nTimeStep; t++ )
		{
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
							nextRange = searchRangeByID( rangeIdList[iR], t );
							assert( nextRange != NULL );
							cout << "here2" << endl;

						}
						else
						{
							//nextRange = rangeBlockTree[t].searchTreeByRangeID( rangeBlockTree[t].getRoot(), rangeIdList[iR] );
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

		}// end for : t

	}// End function

	FEL_range<T>*
	searchRangeByID( int rangeID, int t = 0 )
	{
		typename list<FEL_range<T> >::iterator rangeIter = rangeBlockList.begin();
		for( int iR=0; iR<rangeID; iR++ )
			rangeIter ++;

		// Is there a better/safer way to return pointers?
		return &(*rangeIter);

	}// End function


	void
	estimateDistributionInARegion( int* queryBounds, double* estDist, int t = 0 )
	{
		assert( estDist!= NULL );
		memset( estDist, 0, sizeof(double)*nBin );
		int queryBoundsLocal[6];
		double estDistLocal[nBin];
		memset( estDistLocal, 0, sizeof(double)*nBin );

		if( encodingType == REGULAR )
		{
			memcpy( queryBoundsLocal, queryBounds, sizeof(int)*6 );
			estimateDistributionInARegion_Regular( queryBoundsLocal, estDistLocal, t );
			memcpy( estDist, estDistLocal, sizeof(double)*nBin );

		} // end if : regular encoding
		else
		{
			memcpy( queryBoundsLocal, queryBounds, sizeof(int)*6 );
			estimateDistributionInARegion_Treebased( queryBoundsLocal, estDistLocal, t );
			//#ifdef DEBUG_MODE
			double m = ITL_util<double>::Min( estDist, nBin );
			double M = ITL_util<double>::Max( estDist, nBin );
			cout << "m1: " << m << " " << "M1: " << M << endl;
			//#endif
			memcpy( estDist, estDistLocal, sizeof(double)*nBin );

		}// end else : tree based encoding

		// Normalize estimated distribution
		double sum = 0;
		for( int i=0; i<nBin; i++ )
			sum += estDist[i];
		for( int i=0; i<nBin; i++ )
			estDist[i] = estDist[i] / sum;

	}// End function

	void
	estimateDistributionInARegion_Regular( int* queryBounds, double* estDist, int t = 0 )
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

	}// End function

	void
	estimateDistributionInARegion_Treebased( int* queryBounds, double* estDist, int t = 0 )
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
		overlappingLeafRangeList = rangeBlockTree[t].searchLeafNodesByRegion( rangeBlockTree[t].getRoot(), queryBoundsLocal );
		fprintf( stderr, "%d leaf nodes overlap with the query ...\n", overlappingLeafRangeList.size() );

		// Scan through the list
		// of returned ranges and collect distributions
		//for( typename list<FEL_range<T>* >::iterator rangeIter = overlappingLeafRangeList.begin();
		for( typename list<FEL_range<T> >::iterator rangeIter = overlappingLeafRangeList.begin();
			 rangeIter != overlappingLeafRangeList.end();
			 ++rangeIter )
		{
			// Get next leaf range
			//leafRange = *rangeIter;

			// Get range bounds
			rangeIter->getBounds( rangeBounds );
			//leafRange->getBounds( rangeBounds );

			// Get ID of matching domain
			matchedDomainID = rangeIter->getMatchedDomainId();
			//matchedDomainID = leafRange->getMatchedDomainId();

			#ifdef DUBUG_MODE
			fprintf( stderr, "Query [<%d,%d,%d>-<%d,%d,%d>] overlaps with range [<%d,%d,%d>-<%d,%d,%d>]\n",
					 queryBounds[0], queryBounds[1], queryBounds[2],
					 queryBounds[3], queryBounds[4], queryBounds[5],
					 rangeBounds[0], rangeBounds[1], rangeBounds[2],
					 rangeBounds[3], rangeBounds[4], rangeBounds[5] );
			fprintf( stderr, "Matched Domain: %d\n", matchedDomainID );
			#endif

			// Get distribution and transformation of matching domain
			typename list<FEL_domain<T> >::iterator domIter = domainList.begin();
			for( int iD=0; iD<matchedDomainID; iD++ )
				domIter++;
			domIter->getDistribution( matchedDomainFreqList, &nb );
			rangeIter->getMatchTransformation( nRot, isRef );
			//leafRange->getMatchTransformation( nRot, isRef );
			//cout << nRot << " " << isRef << endl;

			// Transform domain distribution to match range
			transformDistribution( matchedDomainFreqList, nBin, nRot, isRef );

			// Update estimated distribution
			ITL_util<double>::addArrays( estDistLocal, matchedDomainFreqList, estDistLocal, nBin );
			#ifdef DEBUG_MODE
			double m = ITL_util<double>::Min( estDistLocal, nBin );
			double M = ITL_util<double>::Max( estDistLocal, nBin );
			cout << "m1: " << m << " " << "M1: " << M << endl;
			#endif

		}// end for : overlappingLeafRangeList

		memcpy( estDist, estDistLocal, sizeof(double)*nBin );

	}// End function


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

		double buffer;
		double buffer2[nBin-1];
		if( nRot > 0 )
		{
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
		}

		memcpy( freqList, freqListTemp, sizeof(double)*nBin );

		// Free temporary resources
		delete [] freqListTemp;
	}


	int
	getNumTotalDomain()
	{
		return nTotalDomain;

	}// End function

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
	}

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
	}

	void
	getDomainBlockTimevaryingEntropy( int id, float *entropyList )
	{
		assert( entropyList != NULL );

		typename list<FEL_domain<T> >::iterator tmpIter = domainList.begin();
		for( int iD=0; iD<id; iD++ )
			tmpIter++;

		float localEntropyList[nTimeStep];
		tmpIter->getTimevaryingEntropy( localEntropyList );
		memcpy( entropyList, localEntropyList, sizeof(float)*nTimeStep );
	}

	void
	getDomainBlockTimevaryingDistribution( int id, double *fList )
	{
		assert( fList != NULL );

		typename list<FEL_domain<T> >::iterator tmpIter = domainList.begin();
		for( int iD=0; iD<id; iD++ )
			tmpIter++;

		double tmp[nBin*nTimeStep];
		int nb, nt;
		tmpIter->getTimeVaryingDistribution( tmp, &nb, &nt );
		memcpy( fList, tmp, sizeof(double)*(nBin*nTimeStep) );

		#ifdef DEBUG_MODE
		double m = ITL_util<double>::Min( tmp, nBin*nTimeStep );
		double M = ITL_util<double>::Max( tmp, nBin*nTimeStep );
		cout << "m1: " << m << " " << "M1: " << M << endl;
		#endif


	}

	void
	getDomainBlockTimevaryingDistribution_Part( int* domainid, int* segid, double *fList )
	{
		assert( fList != NULL );
		int tstart = 0;
		int ntPerSeg = nTimeStep / nSegment;
		double tmp[nBin*nTimeStep];

		for( int iS=0; iS<nSegment; iS++ )
		{
			typename list<FEL_domain<T> >::iterator tmpIter = domainList.begin();
			for( int iD=0; iD<domainid[iS]; iD++ )
				tmpIter++;

			int nb, nt;
			tmpIter->getTimeVaryingDistributionOfSegment( iS*ntPerSeg, ntPerSeg, tmp+nBin*(iS*ntPerSeg), &nb, &nt );
		}

		memcpy( fList, tmp, sizeof(double)*(nBin*nTimeStep) );

		#ifdef DEBUG_MODE
		double m = ITL_util<double>::Min( tmp, nBin*nTimeStep );
		double M = ITL_util<double>::Max( tmp, nBin*nTimeStep );
		cout << "m1: " << m << " " << "M1: " << M << endl;
		#endif


	}

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

	FEL_spacetree<T>*
	getSpaceTree( int t )
	{
		return ( &rangeBlockTree[t] );
	}

	int getNumTotalRange()
	{
		return nTotalRange;
	}

	//FEL_range<VECTOR3> getRangeBlock( int id );
	void
	getRangeBlockTimevaryingDistribution( int rangeid, double* fList, int* matchingDomainID,
										  double* error, int* optRot, bool* optRef )
	{
		assert( fList != NULL );
		double tmp[nBin*nTimeStep];
		int nb, nt;

		if( encodingType == REGULAR )
		{
			fprintf( stderr, "Searching range block list for the range block with ID: %d\n", rangeid );
			// Linearly scan the range list to find the queried one
			// if regular encoding was used
			typename list<FEL_range<T> >::iterator tmpIter = rangeBlockList.begin();
			for( int iD=0; iD<rangeid; iD++ )
				tmpIter++;

			// Get distribitution
			tmpIter->getTimevaryingDistribution( tmp, &nb, &nt );

			// Get matching domain ID
			(*matchingDomainID) = tmpIter->getMatchedDomainId();

			// Get Match error
			(*error) = tmpIter->getMatchError();

			// Get optimal transformation
			tmpIter->getMatchTransformation( *optRot, *optRef );

			// Copy
			memcpy( fList, tmp, sizeof(double)*(nBin*nTimeStep) );
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
			rangePtr->getTimevaryingDistribution( tmp, &nb, &nt );

			//cout << "2" << endl;

			// Get matching domain
			(*matchingDomainID) = rangePtr->getMatchedDomainId();

			// Get Match error
			(*error) = rangePtr->getMatchError();

			// Get optimal transformation
			rangePtr->getMatchTransformation( *optRot, *optRef );

			//cout << "3" << endl;

			// Copy
			memcpy( fList, tmp, sizeof(double)*(nBin*nTimeStep) );

		}

		//double m = ITL_util<double>::Min( tmp, nBin );
		//double M = ITL_util<double>::Max( tmp, nBin );
		//cout << "m1: " << m << " " << "M1: " << M << endl;
	}

	//FEL_range<VECTOR3> getRangeBlock( int id );
	void
	getRangeBlockTimevaryingDistribution_Part( int rangeid,
											   double* fList,
											   int* matchingDomainID,
											   int* matchingSegID,
										  	   double* error, int* optRot, bool* optRef )
	{
		assert( fList != NULL );
		double tmp[nBin*nTimeStep];
		int matchedDomainIdList[nSegment];
		int matchedSegmentIdList[nSegment];
		double matchErrorList[nSegment];
		int optRotList[nSegment];
		bool optRefList[nSegment];

		int nb, nt;

		if( encodingType == REGULAR )
		{
			fprintf( stderr, "Searching range block list for the range block with ID: %d\n", rangeid );
			// Linearly scan the range list to find the queried one
			// if regular encoding was used
			typename list<FEL_range<T> >::iterator tmpIter = rangeBlockList.begin();
			for( int iD=0; iD<rangeid; iD++ )
				tmpIter++;

			// Get time varying distribitution
			tmpIter->getTimevaryingDistribution( tmp, &nb, &nt );

			// Get matching domain ID
			tmpIter->getMatchedDomainId( matchedDomainIdList, matchedSegmentIdList );

			// Get Match error
			tmpIter->getMatchError( matchErrorList );

			// Get optimal transformation
			tmpIter->getMatchTransformation( optRotList, optRefList );

		}
		else
		{
			fprintf( stderr, "Searching range block tree for the range block with ID: %d\n", rangeid );
			fprintf( stderr, "**********Not implemented yet**************\n" );

			/*
			// Search quadtree/octree to find the queried one
			// if hierarchical encoding was used
			FEL_range<T>* rangePtr = rangeBlockTree.searchTreeByRangeID( rangeBlockTree.getRoot(), rangeid );

			//cout << "1" << endl;
			assert( rangePtr != NULL );

			// Get distribitution
			rangePtr->getTimevaryingDistribution( tmp, &nb, &nt );

			//cout << "2" << endl;

			// Get matching domain
			(*matchingDomainID) = rangePtr->getMatchedDomainId();

			// Get Match error
			(*error) = rangePtr->getMatchError();

			// Get optimal transformation
			rangePtr->getMatchTransformation( *optRot, *optRef );

			//cout << "3" << endl;

			// Copy
			memcpy( fList, tmp, sizeof(double)*(nBin*nTimeStep) );
			*/

		}

		// Copy
		memcpy( fList, tmp, sizeof(double)*(nBin*nTimeStep) );
		memcpy( matchingDomainID, matchedDomainIdList, sizeof(int)*nSegment );
		memcpy( matchingSegID, matchedSegmentIdList, sizeof(int)*nSegment );
		memcpy( error, matchErrorList, sizeof(double)*nSegment );
		memcpy( optRot, optRotList, sizeof(int)*nSegment );
		memcpy( optRef, optRefList, sizeof(bool)*nSegment );

		//double m = ITL_util<double>::Min( tmp, nBin );
		//double M = ITL_util<double>::Max( tmp, nBin );
		//cout << "m1: " << m << " " << "M1: " << M << endl;
	}

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

	}

	void filterRangeBlockExtents( int lod, int* nReturnedRanges, float* rangeExtentList, int t = 0 )
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
				rangeAtLodList = rangeBlockTree[t].searchTreeByLevel( rangeBlockTree[t].getRoot(), lod );
			else if( lod == -1 )
				rangeAtLodList = rangeBlockTree[t].getLeafNodes( rangeBlockTree[t].getRoot() );

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
	}

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

				matlab_plotting_on = false;
			}
		}
	}// end function

	void
	plotEncodingResutls()
	{
		// Draw Matlab plots
		if( matlab_plotting_on )
		{
			plotter->displayDomainUtilization_Matlab();
			plotter->displayMatchingPairEntropyLevel_Matlab();
			plotter->displayMatchErrorDistribution_Matlab();

		}
	}

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

	void printRanges( int t )
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
		}
	}

	~FEL_encoder_timevarying()
	{
		domainList.clear();

		//if( inputBinField != NULL )	delete [] inputBinField;
		//if( downsampledDataField != NULL ) delete [] downsampledDataField;
		//if( downsampledBinField != NULL ) delete [] downsampledBinField;

		//for( int t=0; t<nTimeStep; t++ )
		//{
		//	if( downsampledDataFieldArray[t] != NULL ) delete [] downsampledDataFieldArray[t];
		//	if( downsampledBinFieldArray[t] != NULL ) delete [] downsampledBinFieldArray[t];
		//}

		//if( downsampledDataFieldArray != NULL )	delete [] downsampledDataFieldArray;
		//if( downsampledBinFieldArray != NULL )	delete [] downsampledBinFieldArray;

	}
};

#endif
/* FEL_ENCODER_H_ */
