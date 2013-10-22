/*
 * FEL_encoder_incremental.cpp
 *
 *  Created on: Sep 10, 2012
 *      Author: abon
 */

#include "FEL_encoder_incremental.h"

FEL_encoder_incremental::FEL_encoder_incremental()
{
	nBin = 0;
	nDomain = 0;
	nUsedDomain = 0;
	minDomainEntropy = maxDomainEntropy = 0;
	distanceMeasureType = 1;	// L2-error
	matchErrorThreshold = 0;

	memset( domainUtilizationList, 0, sizeof(int)*50000 );

}// end constructor


FEL_encoder_incremental::FEL_encoder_incremental( int nbin, list<FEL_domain_core>* domainlist,
														 double threshold )
{
	nBin = nbin;
	nDomain = domainlist->size();
	nUsedDomain = 0;
	minDomainEntropy = maxDomainEntropy = 0;
	distanceMeasureType = 1;	// L2-error
	matchErrorThreshold = threshold;

	memset( domainUtilizationList, 0, sizeof(int)*50000 );

	list<FEL_domain_core>::iterator iter = domainList.begin();
	domainList.insert( iter, domainlist->begin(), domainlist->end() );

}// end constructor

FEL_encoder_incremental::FEL_encoder_incremental( const FEL_encoder_incremental& that )
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

		memcpy( this->domainUtilizationList, that.domainUtilizationList, sizeof(int)*50000 );
	}

}// end copy constructor

FEL_encoder_incremental&
FEL_encoder_incremental::operator= ( const FEL_encoder_incremental& that )
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

		memcpy( this->domainUtilizationList, that.domainUtilizationList, sizeof(int)*50000 );

	}
	// by convention, always return *this
	return *this;
}// end assignment operator

FEL_encoder_incremental::~FEL_encoder_incremental()
{
	if( domainList.size() > 0 )	domainList.clear();
}// end destructor

void
FEL_encoder_incremental::init( int nbin, list<FEL_domain_core>* domainlist,
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

int
FEL_encoder_incremental::encode_Distribution( double* rangeFreqList, long nDataPoint,
					 	 	 	 	 	 	  FEL_encodedhistogram* encodedHist,
					 	 	 	 	 	 	  double errorThreshold,
					 	 	 	 	 	 	  bool allZeroFlag )
{
	int optRot = 0, bestMatch = -1, nHighErrorBin = 0;
	bool optRef = false;
	double matchError = -1;
	int zeroDomainId = nDomain-1;
	double warningThreshold = 0.75;

	#ifdef DEBUG_MODE
	double m = ITL_util<double>::Min( rangeFreqList, nBin );
	double M = ITL_util<double>::Max( rangeFreqList, nBin );
	fprintf( stderr, "in encode: Range distribution range: %g, %g\n", m, M );
	#endif

	// Handle special case - all zero distribution
	if( allZeroFlag == true )
	{
		bestMatch = zeroDomainId;
		matchError = 0;
		optRot = 0;
		optRef = false;

		// Keep track of utilized domains
		incrementDomainUtilization( bestMatch );
	}
	else
	{
		// Iterate though all the domains
		// and find the nearest one
		bestMatch = matchRangeToDomains2( rangeFreqList, &matchError, &optRot, &optRef );
		//bestMatch = matchRangeToDomains_Entropybased( rangeFreqList, &matchError, &optRot, &optRef );

		//if( bestMatch == -1 )
		//{
		//	double limits[6] = { 0, 0, 0, 1, 1, 1 };
		//	FEL_domain_core newDomain( limits, nBin, rangeFreqList );
		//	domainList.push_back( newDomain );
		//	domainTree.
		//	nDomain ++;

		//	bestMatch = nDomain-1;
		//	matchError = 0;
		//	optRot = 0;
		//	optRef = false;
		//}

		// Keep track of utilized domains
		incrementDomainUtilization( bestMatch );
	}

	//fprintf( stderr, "from active 1: Total active count: %d\n", nUsedDomain );

	#ifdef DEBUG_MODE
	fprintf( stderr, "Range maps to %d-th domain with %g <%d, %d> error ...\n",
			 	     bestMatch, matchError, optRot, (int)optRef );
	#endif

	// Encode range distribution
	encodedHist->setMatchingDomainId( bestMatch );
	encodedHist->setNumDataElement( nDataPoint );
	encodedHist->setRotation( optRot );
	encodedHist->setReflection( optRef );

	//fprintf( stderr, "from active 2: Total active count: %d\n", nUsedDomain );

	if( errorThreshold != 0 && matchError > errorThreshold )
	{
		double* domainDist = new double[nBin];
		double* transformedDomain = new double[nBin];
		int* highErrorBinIds = new int[nBin];
		double* highErrors = new double[nBin];

		// Get domain distribution
		list<FEL_domain_core>::iterator iter = domainList.begin();
		for( int i=0; i< bestMatch; i++ )
			iter++;
		iter->getDistribution( domainDist );

		// Transform domain
		memcpy( transformedDomain, domainDist, sizeof(double)*nBin );
		// ****transformDistribution( transformedDomain, optRot, optRef );

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

		// Store high error bins
		encodedHist->setNumHighErrorBin( nHighErrorBin );
		encodedHist->setHighErrorBinInfo( highErrorBinIds, highErrors );

		delete [] highErrors;
		delete [] highErrorBinIds;
		delete [] transformedDomain;
		delete [] domainDist;

	}// end if

	//fprintf( stderr, "from active 3: Total active count: %d\n", nUsedDomain );

	if( matchError > errorThreshold &&
		nHighErrorBin > warningThreshold*nBin )
	{
		return 1;
	}
	return 0;
}

int
FEL_encoder_incremental::encode_Distribution_2( double* spanFreqList,
													  double* residualFreqList,
													  long nDataPoint,
					 	 	 	 	 	 	  	  	  FEL_encodedhistogram* encodedHist )
{
	int optRot = 0, bestMatch = -1;
	bool optRef = false;
	double residue[nBin];

	#ifdef DEBUG_MODE
	double m = ITL_util<double>::Min( spanFreqList, nBin );
	double M = ITL_util<double>::Max( spanFreqList, nBin );
	fprintf( stderr, "in encode: integral distribution range: %g, %g\n", m, M );
	#endif

	// Iterate though all the templates
	// and find the nearest one
	bestMatch = findBestMatch_peakbased( spanFreqList,
										 residue,
			 	 	 	 	 	 	 	 &optRot,
			 	 	 	 	 	 	 	 &optRef );
	memcpy( residualFreqList, residue, sizeof(double)*nBin );

	#ifdef DEBUG_MODE
	fprintf( stderr, "Range maps to %d-th domain <%d, %d> error ...\n",
			 	     bestMatch, optRot, (int)optRef );
	#endif

	// Encode range distribution
	encodedHist->setMatchingDomainId( bestMatch );
	encodedHist->setNumDataElement( nDataPoint );
	encodedHist->setRotation( optRot );
	encodedHist->setReflection( optRef );

	//fprintf( stderr, "from active 2: Total active count: %d\n", nUsedDomain );
}

int
FEL_encoder_incremental::encode_Distribution_2( double* spanFreqList,
													  double* residualFreqList,
													  long nDataPoint,
					 	 	 	 	 	 	  	  	  int* optTemplateId,
					 	 	 	 	 	 	  	  	  int* optRotation,
					 	 	 	 	 	 	  	  	  bool* reflectionFlag )
{
	int optRot = 0, bestMatch = -1;
	bool optRef = false;
	double residue[nBin];

	#ifdef DEBUG_MODE
	double m = ITL_util<double>::Min( spanFreqList, nBin );
	double M = ITL_util<double>::Max( spanFreqList, nBin );
	fprintf( stderr, "in encode: integral distribution range: %g, %g\n", m, M );
	#endif

	// Iterate though all the templates
	// and find the nearest one
	bestMatch = findBestMatch_peakbased( spanFreqList,
										 residue,
			 	 	 	 	 	 	 	 &optRot,
			 	 	 	 	 	 	 	 &optRef );
	memcpy( residualFreqList, residue, sizeof(double)*nBin );

	//#ifdef DEBUG_MODE
	//fprintf( stderr, "Range maps to %d-th domain <%d, %d>\n",
	//		 	     bestMatch, optRot, (int)optRef );
	//#endif

	// Encode range distribution
	(*optTemplateId) = bestMatch;
	(*optRotation) = optRot;
	(*reflectionFlag) = optRef;

}


void
FEL_encoder_incremental::correctBadEncoding(  FEL_encodedhistogram* encodedHist,
													double* rangeFreqList, int nData )
{
	double newDomainLimits[6] = {0, 1, 0, 1, 0, 1};

	FEL_domain_core newDomain( newDomainLimits, nBin, rangeFreqList );
	addNewDomain( newDomain );

	// Decrement the utilization count of current domain
	decrementDomainUtilization( encodedHist->getMatchingDomainId() );

	// Update encoding of curent range distribution
	encodedHist->setMatchingDomainId( nDomain-1 );
	encodedHist->setNumDataElement( nData );
	encodedHist->setRotation( 0 );
	encodedHist->setReflection( false );
	encodedHist->setNumHighErrorBin( 0 );

	// Decrement the utilization count of new domain
	incrementDomainUtilization( encodedHist->getMatchingDomainId() );
}

int
FEL_encoder_incremental::matchRangeToDomains2( double* rangeFreqList, double *error,
											   	   	 int* rotateAmount, bool* isReflected )
{
	int domainId = 0, bestMatchedDomainId = -1;
	double e = 0, minError = 1000000;
	double domainFreqList[nBin];
	double bestDomainFreqList[nBin];
	int rotAmount = 0;
	bool isRef = false;

	// Scan through all domain distirbutions
	for( list<FEL_domain_core>::iterator domainIter = domainList.begin();
		 domainIter != domainList.end();
		 domainIter++ )
	{


		// Compare only if the domain is active
		if( domainIter->getActiveStatus() )
		{

			// Get distribution of next domain
			domainIter->getDistribution( domainFreqList );

			// Compute match error between range and domain distributions
			e =  FEL_util<float>::matchRangeToDomainByDistribution( distanceMeasureType,
																	rangeFreqList, domainFreqList,
																	nBin,
																	&rotAmount, &isRef );

			// Initialize min error
			//if( domainId == 0 )	minError = e;

			//cout << "domainID: " << domainId << " e: " << e << endl;
			if( e < minError )
			{
				minError = e;
				(*error) = minError;
				(*rotateAmount) = rotAmount;
				(*isReflected) = isRef;
				bestMatchedDomainId = domainId;
				memcpy( bestDomainFreqList, domainFreqList, sizeof( double )*nBin );
			}
		}

		// Move to next domain
		domainId ++;
	}

	return bestMatchedDomainId;
}

int
FEL_encoder_incremental::findBestMatch_peakbased( double* freqList,
														 double* errorList,
														 int* rotateAmount,
														 bool* isReflected )
{
	int domainId = 0, bestMatchedDomainId = -1;
	double e = 0, minError = 1000000;
	double templateFreqList[nBin];
	double bestDomainFreqList[nBin];
	double errorFreqList[nBin];
	int rotAmount = 0;
	bool isRef = false;

	// Scan through all templates
	for( list<FEL_domain_core>::iterator domainIter = domainList.begin();
		 domainIter != domainList.end();
		 domainIter++ )
	{
		// Get template distribution
		domainIter->getDistribution( templateFreqList );

		// Compute match error between range and domain distributions
		FEL_util<double>::compareDistributionWithTemplate( freqList, templateFreqList,
															nBin,
							 	 	 	 	 	 	 	    &e, errorFreqList,
							 	 	 	 	 	 	 	    &rotAmount, &isRef );

		if( e < minError )
		{
			minError = e;
			(*rotateAmount) = rotAmount;
			(*isReflected) = isRef;
			bestMatchedDomainId = domainId;
			memcpy( errorList, errorFreqList, sizeof(double)*nBin );
		}

		// Move to next domain
		domainId ++;

	}// end for: scan through templated

	return bestMatchedDomainId;
}

int
FEL_encoder_incremental::matchRangeToDomains_Entropybased( double* rangeFreqList, double *error,
														   int* rotateAmount, bool* isReflected  )
{
	int domainId = 0, bestMatchedDomainId = -1;
	double e, minError;
	double domainFreqList[nBin];
	double bestDomainFreqList[nBin];
	int rotAmount = 0;
	bool isRef = false;
	double rangeEntropy;
	double delH, entropyRangeLow, entropyRangeHigh;
	list<int> searchResultDomainList;

	// Compute entropy of the query distribution
	rangeEntropy = ITL_entropycore::computeEntropy_HistogramBased2( rangeFreqList, nBin, false );
	#ifdef DEBUG_MODE
	fprintf( stderr, "Entropy of range distribution: %g\n", rangeEntropy );
	#endif

	// Create a query entropy range
	// Keep expanding the query unless
	// a domain is found in the range
	double maxEntropy = ( log(nBin) / log(2) );
	double percnt = 0.1;
	do
	{
		delH = ( percnt * maxEntropy );
		entropyRangeLow = rangeEntropy - delH;
		entropyRangeHigh = rangeEntropy + delH;
		#ifdef DEBUG_MODE
		fprintf( stderr, "Entropy range for querying: <%g, %g>\n", entropyRangeLow, entropyRangeHigh );
		#endif

		// Find domains within query range
		searchResultDomainList = domainTree.searchTree( domainTree.getRoot(), entropyRangeLow, entropyRangeHigh );

		// Update range
		percnt = percnt + 0.05;

	}while( searchResultDomainList.size() == 0 && percnt <= 1 );

	if( searchResultDomainList.size() == 0 )
	{
		fprintf( stderr, "This shouldn't happen !!\n" );
		exit(0);
	}

	// Match range with all returned domains
	for( list<int>::iterator it = searchResultDomainList.begin();
		 it != searchResultDomainList.end();
		 it++ )
	{
		// Retrieve next domain within entropy range
		int nextDomainId = (*it);
		list<FEL_domain_core>::iterator domainIter = domainList.begin();
		for( int iD=0; iD<nextDomainId; iD++ )
			domainIter++;

		// Get distribution of next domain
		domainIter->getDistribution( domainFreqList );

		// Compute match error between range and domain distributions
		e =  FEL_util<double>::matchRangeToDomainByDistribution( distanceMeasureType,
																 rangeFreqList, domainFreqList,
																 nBin,
																 &rotAmount, &isRef );

		// Initialize min error
		if( domainId == 0 )	minError = e;

		if( e < minError )
		{
			minError = e;
			(*error) = minError;
			(*rotateAmount) = rotAmount;
			(*isReflected) = isRef;
			bestMatchedDomainId = domainId;
			memcpy( bestDomainFreqList, domainFreqList, sizeof( double )*nBin );
		}

		domainId ++;
	}

	return bestMatchedDomainId;
}

void
FEL_encoder_incremental::updateDomainSet()
{
	int domainId = 0;

	for( list<FEL_domain_core>::iterator domainIter = domainList.begin();
		 domainIter != domainList.end();
		 domainIter++ )
	{
		if( domainUtilizationList[domainId] == 0 )
			domainIter->setActiveStatus( false );

		domainId ++;
	}

	#ifdef DEBUG_MODE
	fprintf( stderr, "Total domain count: %d\n", nDomain );
	fprintf( stderr, "Active domain count: %d\n", nUsedDomain );
	#endif
}

void
FEL_encoder_incremental::addNewDomain( FEL_domain_core newDomain )
{
	newDomain.setActiveStatus( true );
	domainList.push_back( newDomain );
	nDomain = nDomain + 1;

	#ifdef DEBUG_MODE
	fprintf( stderr, "Total domain count: %d\n", nDomain );
	#endif
}

void
FEL_encoder_incremental::computeDomainEntropies()
{
	int domainBounds[6];
	int nPointDomain = 1;
	double domainFreqList[nBin];
	float domainEntropy;
	int nb;

	int domainID = 0;
	for( list< FEL_domain_core>::iterator domainIter = domainList.begin();
		 domainIter != domainList.end();
		 domainIter++ )
	{
		// Get domain distribution
		domainIter->getDistribution( domainFreqList );

		// Compute and store domain entropy
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

		// Move to next domain
		domainID ++;
	}
}

void
FEL_encoder_incremental::sortDomains_Entropybased()
{
	// Sort domains based on entropy
	domainList.sort( FEL_domain_core::compare_entropy );
}

/*
void
FEL_encoder_incremental::pruneDomains_Entropybased()
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
FEL_encoder_incremental::pruneDomains_Entropybased2()
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
*/

void
FEL_encoder_incremental::organizeDomains( bool isDomainPruningOn )
{
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
		//pruneDomains_Entropybased2();
		fprintf( stderr, "Number of domains after pruning: %d\n", domainList.size() );

		// Update the total number of domains
		nDomain = domainList.size();
	}

	// Organize sorted domains in a segment tree
	fprintf( stderr, "%s: %d: Storing domains in a hierarchy ...\n", __FILE__, __LINE__ );
	domainTree.initTree( domainList );
	domainTree.createTree( domainTree.getRoot(), 0, nDomain-1, 0 );
}

FEL_segmenttree2*
FEL_encoder_incremental::getDomainTree()
{
	return (&domainTree);
}

int
FEL_encoder_incremental::getDomainCount()
{
	return nDomain;
}

int
FEL_encoder_incremental::getActiveDomainCount()
{
	return nUsedDomain;
}

void
FEL_encoder_incremental::incrementDomainUtilization( int id )
{
	//assert( id < 50000 );
	if( domainUtilizationList[id] == 0 )
		nUsedDomain ++;
	domainUtilizationList[id] ++;
}

void
FEL_encoder_incremental::decrementDomainUtilization( int id )
{
	assert( id < 10000 );
	domainUtilizationList[id] --;
	if( domainUtilizationList[id] == 0 )
		nUsedDomain --;
	if( domainUtilizationList[id] < 0 )
		domainUtilizationList[id] = 0;
}

