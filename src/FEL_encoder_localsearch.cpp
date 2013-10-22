/*
 * FEL_encoder_localsearch.cpp
 *
 *  Created on: Mar 20, 2013
 *      Author: abon
 */

#include "FEL_encoder_localsearch.h"

FEL_encoder_localsearch::FEL_encoder_localsearch()
{
	nBin = 0;
	nDomain = 0;
	nUsedDomain = 0;
	minDomainEntropy = maxDomainEntropy = 0;
	distanceMeasureType = 1;	// L2-error
	matchErrorThreshold = 0;

}// end constructor


FEL_encoder_localsearch::FEL_encoder_localsearch( int nbin, list<FEL_domain_core>* domainlist,
														 double threshold )
{
	nBin = nbin;
	nDomain = domainlist->size();
	nUsedDomain = 0;
	minDomainEntropy = maxDomainEntropy = 0;
	distanceMeasureType = 1;	// L2-error
	matchErrorThreshold = threshold;

	list<FEL_domain_core>::iterator iter = domainList.begin();
	domainList.insert( iter, domainlist->begin(), domainlist->end() );

}// end constructor

FEL_encoder_localsearch::FEL_encoder_localsearch( const FEL_encoder_localsearch& that )
{
	this->nBin = that.nBin;
	this->nDomain = that.nDomain;
	this->nUsedDomain = that.nUsedDomain;
	this->minDomainEntropy = that.minDomainEntropy;
	this->maxDomainEntropy = that.maxDomainEntropy;
	this->distanceMeasureType = that.distanceMeasureType;
	this->matchErrorThreshold = that.matchErrorThreshold;

	if( that.domainList.size() > 0 )
	{
		list<FEL_domain_core>::iterator iter = this->domainList.begin();
		this->domainList.insert( iter, that.domainList.begin(), that.domainList.end() );
	}

}// end copy constructor

FEL_encoder_localsearch&
FEL_encoder_localsearch::operator= ( const FEL_encoder_localsearch& that )
{
	if (this != &that ) // protect against invalid self-assignment
	{
		this->nBin = that.nBin;
		this->nDomain = that.nDomain;
		this->nUsedDomain = that.nUsedDomain;
		this->minDomainEntropy = that.minDomainEntropy;
		this->maxDomainEntropy = that.maxDomainEntropy;
		this->distanceMeasureType = that.distanceMeasureType;
		this->matchErrorThreshold = that.matchErrorThreshold;

		if( that.domainList.size() > 0 )
		{
			list<FEL_domain_core>::iterator iter = this->domainList.begin();
			this->domainList.insert( iter, that.domainList.begin(), that.domainList.end() );
		}
	}
	// by convention, always return *this
	return *this;
}// end assignment operator

FEL_encoder_localsearch::~FEL_encoder_localsearch()
{
	if( domainList.size() > 0 )	domainList.clear();
}// end destructor

void
FEL_encoder_localsearch::init( int nbin, list<FEL_domain_core>* domainlist,
								   double threshold )
{
	nBin = nbin;
	nDomain = domainlist->size();
	nUsedDomain = 0;
	minDomainEntropy = maxDomainEntropy = 0;
	distanceMeasureType = 0;
	matchErrorThreshold = threshold;

	list<FEL_domain_core>::iterator iter = domainList.begin();
	domainList.insert( iter, domainlist->begin(), domainlist->end() );

}// end function

void
FEL_encoder_localsearch::organizeDomains()
{
	double domainCenter[3];

	// Create data structure for storing reduced domains
	domainCenterList = vtkSmartPointer<vtkPoints>::New();

	// Iterate over all domains
	vtkIdType pId = 0;
	for( list<FEL_domain_core>::iterator iter = domainList.begin();
		 iter != domainList.end();
		 iter++ )
	{
		// Compute domain center
		iter->getCenter( domainCenter );

		// Store coefficients in quadtree (octree)
		domainCenterList->InsertPoint( pId, domainCenter );

		// Move to next domain
		pId ++;

	}// end for

	// Create polydata
	fprintf( stderr, "Creating polydata ...\n" );
	domainPolyData = vtkSmartPointer<vtkPolyData>::New();
	domainPolyData->SetPoints( domainCenterList );

	//Create the octree
	fprintf( stderr, "Creating octree ...\n" );
	octree = vtkSmartPointer<vtkOctreePointLocator>::New();
	octree->SetDataSet( domainPolyData );
	octree->BuildLocator();

}


void
FEL_encoder_localsearch::encode_Distribution_LocalSearch( double* rangeFreqList, double* rangeLimits,
		 	 	 	 	 	 	 	 	 	 	 	 	 	 	  double distanceThreshold,	long nDataPoint,
					 	 	 	 	 	 	  	  	  	  	  	  FEL_encodedhistogram* encodedHist,
					 	 	 	 	 	 	  	  	  	  	  	  double errorThreshold, bool allZeroFlag )
{
	int optRot = 0, bestMatch = -1, nHighErrorBin = 0;
	int lowSub[3], highSub[3];
	bool optRef = false;
	double matchError = -1;
	double rangeCenter[3];
	double domainDist[nBin];
	double transformedDomain[nBin];
	int highErrorBinIds[nBin];
	double highErrors[nBin];

	#ifdef DEBUG_MODE
	double m = ITL_util<double>::Min( rangeFreqList, nBin );
	double M = ITL_util<double>::Max( rangeFreqList, nBin );
	fprintf( stderr, "in encode: Range distribution range: %g, %g\n", m, M );
	#endif

	// Compute range center
	// Get limit and center of span
	lowSub[0] = rangeLimits[0]; highSub[0] = rangeLimits[1];
	lowSub[1] = rangeLimits[2]; highSub[1] = rangeLimits[3];
	lowSub[2] = rangeLimits[4]; highSub[2] = rangeLimits[5];
	rangeCenter[0] = ( rangeLimits[0] + rangeLimits[1] ) / 2.0f;
	rangeCenter[1] = ( rangeLimits[2] + rangeLimits[3] ) / 2.0f;
	rangeCenter[2] = ( rangeLimits[4] + rangeLimits[5] ) / 2.0f;

	// Iterate though the spatially nearby domains
	bestMatch = matchRangeToDomains_LocalSearch( rangeFreqList, rangeCenter,
												 distanceThreshold,
												 &matchError, &optRot, &optRef );
	#ifdef DEBUG_MODE
	fprintf( stderr, "Current Threshold: %g, best match: %d\n", distanceThreshold, bestMatch );
	fprintf( stderr, "Range maps to %d-th domain with %g <%d, %d> error ...\n",
			 	     bestMatch, matchError, optRot, (int)optRef );
	#endif

	// Case 1: The range did not find any match
	// Add the current range to the domain pool
	//if( bestMatch == -1 )
	//{
	//	bestMatch = correctBadEncoding( rangeFreqList, rangeLimits, &matchError, &optRot, &optRef );
		#ifdef DEBUG_MODE
		fprintf( stderr, "After correction: Range maps to %d-th domain with %g <%d, %d> error ...\n",
				 	     	 	 	 	 	 	 	 	 bestMatch, matchError, optRot, (int)optRef );
		#endif

	//}

	// Case 2: The range did find a match, but with poor accuracy
	// Add the high error bins to reduce error
	if( errorThreshold != 0 && matchError > errorThreshold )
	{
		// Get domain distribution
		list<FEL_domain_core>::iterator iter = domainList.begin();
		for( int i=0; i< bestMatch; i++ )
			iter++;
		iter->getDistribution( domainDist );

		// Transform domain
		memcpy( transformedDomain, domainDist, sizeof(double)*nBin );
		FEL_util<double>::transformDistribution( transformedDomain, nBin, optRot, optRef );

		// Find high error bins and store errors
		int iP = 0;
		for( int i=0; i<nBin; i++ )
		{
			if( abs( transformedDomain[i] - rangeFreqList[i] ) > errorThreshold )
			{
				highErrorBinIds[iP] = i;
				highErrors[iP] = ( rangeFreqList[i] - transformedDomain[i] );

				iP++;
				nHighErrorBin ++;
			}
		}

	}// end if

	// Encode range distribution
	encodedHist->setSpan( lowSub, highSub );
	encodedHist->setMatchingDomainId( bestMatch );
	encodedHist->setNumDataElement( nDataPoint );
	encodedHist->setRotation( optRot );
	encodedHist->setReflection( optRef );
	encodedHist->setEncodingError( matchError );
	if( nHighErrorBin > 0 )
	{
		encodedHist->setNumHighErrorBin( nHighErrorBin );
		encodedHist->setHighErrorBinInfo( highErrorBinIds, highErrors );
	}
}

int
FEL_encoder_localsearch::matchRangeToDomains_LocalSearch( double* rangeFreqList, double* rangeCenter,
																  double distanceThreshold,
																  double *error,
																  int* rotateAmount, bool* isReflected )
{
	int domainId = 0, bestMatchedDomainId = -1;
	double e = 0, minError = 1000000;
	double domainSearchBox[6];
	double domainFreqList[nBin];
	double bestDomainFreqList[nBin];
	int rotAmount = 0;
	bool isRef = false;
	double nextDomainCenter[3];

	#ifdef DEBUG_MODE
	fprintf( stderr, "range Center: %g, %g, %g\n", rangeCenter[0],
													rangeCenter[1],
													rangeCenter[2] );
	#endif
	//for( int ii=0; ii<nBin-1; ii++ )
	//	fprintf( stderr, "%g, ", rangeFreqList[ii] );
	//fprintf( stderr, "%g\n", rangeFreqList[nBin-1] );

	// Create search box around current range
	domainSearchBox[0] = rangeCenter[0] - distanceThreshold;
	domainSearchBox[1] = rangeCenter[0] + distanceThreshold;
	domainSearchBox[2] = rangeCenter[1] - distanceThreshold;
	domainSearchBox[3] = rangeCenter[1] + distanceThreshold;
	domainSearchBox[4] = rangeCenter[2] - distanceThreshold;
	domainSearchBox[5] = rangeCenter[2] + distanceThreshold;

	// Search for domains within this range
	vtkSmartPointer<vtkIdTypeArray> nearbyDomainIdList = vtkSmartPointer<vtkIdTypeArray>::New();
	octree->FindPointsInArea( domainSearchBox, nearbyDomainIdList, true );
	vtkIdType listSize = nearbyDomainIdList->GetNumberOfTuples();
	//fprintf( stderr, "Number of domains near this range: %d\n", (int)listSize );

	// Scan through nearby domains only
	for( vtkIdType i = 0; i < listSize; i++ )
	{
		long domainIndex = (long)nearbyDomainIdList->GetValue( i );
		//fprintf( stderr, "Index of next domain: %d\n", (int)domainIndex );

		// Obtain an iterator to the domain at that position
		std::list<FEL_domain_core>::iterator domainIter = domainList.begin();
		std::advance( domainIter, domainIndex );

		// Get distribution of next domain
		domainIter->getDistribution( domainFreqList );
		//for( int ii=0; ii<nBin-1; ii++ )
		//	fprintf( stderr, "%g, ", domainFreqList[ii] );
		//fprintf( stderr, "%g\n", domainFreqList[nBin-1] );

		// Compute match error between range and domain distributions
		e =  FEL_util<float>::matchRangeToDomainByDistribution( distanceMeasureType,
																rangeFreqList, domainFreqList,
																nBin,
																&rotAmount, &isRef );

		//cout << "domainID: " << domainId << " e: " << e << endl;
		if( e < minError )
		{
			minError = e;
			(*error) = minError;
			(*rotateAmount) = rotAmount;
			(*isReflected) = isRef;
			bestMatchedDomainId = domainIndex;
			memcpy( bestDomainFreqList, domainFreqList, sizeof( double )*nBin );
		}

	}// end loop through nearby domains

	return bestMatchedDomainId;
}

int
FEL_encoder_localsearch::correctBadEncoding(  double* rangeFreqList, double* rangeLimits,
													double *error,
													int* rotateAmount, bool* isReflected )
{
	fprintf( stderr, "Adding this range to domain list ...\n" );

	FEL_domain_core newDomain( rangeLimits, nBin, rangeFreqList );
	addNewDomain( newDomain );

	// Update encoding of information of current range distribution
	(*error) = 0;
	(*rotateAmount) = 0;
	(*isReflected) = false;

	return (nDomain-1);
}

void
FEL_encoder_localsearch::addNewDomain( FEL_domain_core newDomain )
{
	newDomain.setActiveStatus( true );
	domainList.push_back( newDomain );
	nDomain = nDomain + 1;

	#ifdef DEBUG_MODE
	fprintf( stderr, "Total domain count: %d\n", nDomain );
	#endif
}

int
FEL_encoder_localsearch::getDomainCount()
{
	return nDomain;
}
