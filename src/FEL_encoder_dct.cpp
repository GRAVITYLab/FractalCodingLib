/*
 * FEL_encoder_dct.cpp
 *
 *  Created on: March 18, 2013
 *      Author: abon
 */

#include "FEL_encoder_dct.h"

FEL_encoder_dct::FEL_encoder_dct()
{
	nBin = 0;
	nDomain = 0;
	nUsedDomain = 0;
	distanceMeasureType = 1;	// L2-error
	matchErrorThreshold = 0;

	domainBuffer = 0;
	domainBufferRev = 0;

}// end constructor


FEL_encoder_dct::FEL_encoder_dct( int nbin, list<FEL_domain_core>* domainlist,
														 double threshold )
{
	nBin = nbin;
	nDomain = domainlist->size();
	nUsedDomain = 0;
	distanceMeasureType = 1;	// L2-error
	matchErrorThreshold = threshold;

	list<FEL_domain_core>::iterator iter = domainList.begin();
	domainList.insert( iter, domainlist->begin(), domainlist->end() );

	domainBuffer = new double*[(int)domainList.size()];
	domainBufferRev = new double*[(int)domainList.size()];

	iter= domainList.begin();
	double temp[nBin];
	double tempRev[nBin];
	for( int i=0; i<(int)domainList.size(); i++ )
	{
		domainBuffer[i] = new double[nBin*2];
		domainBufferRev[i] = new double[nBin*2];

		iter->getDistribution( temp );

		memcpy( domainBuffer[i], temp, sizeof(double)*nBin );
		memcpy( domainBuffer[i]+nBin, temp, sizeof(double)*nBin );

		for( int i=0; i<nBin; i++ )
			tempRev[i] = temp[nBin-i-1];

		memcpy( domainBufferRev[i], tempRev, sizeof(double)*nBin );
		memcpy( domainBufferRev[i]+nBin, tempRev, sizeof(double)*nBin );

		iter++;
	}

}// end constructor

FEL_encoder_dct::FEL_encoder_dct( const FEL_encoder_dct& that )
{
	this->nBin = that.nBin;
	this->nDomain = that.nDomain;
	this->nUsedDomain = that.nUsedDomain;
	this->distanceMeasureType = that.distanceMeasureType;
	this->matchErrorThreshold = that.matchErrorThreshold;

	if( that.domainList.size() > 0 )
	{
		list<FEL_domain_core>::iterator iter = this->domainList.begin();
		this->domainList.insert( iter, that.domainList.begin(), that.domainList.end() );

		this->domainBuffer = new double*[(int)this->domainList.size()];
		this->domainBufferRev = new double*[(int)this->domainList.size()];

		iter= this->domainList.begin();
		double temp[nBin];
		double tempRev[nBin];
		for( int i=0; i<(int)this->domainList.size(); i++ )
		{
			this->domainBuffer[i] = new double[nBin*2];
			this->domainBufferRev[i] = new double[nBin*2];

			iter->getDistribution( temp );

			memcpy( this->domainBuffer[i], temp, sizeof(double)*nBin );
			memcpy( this->domainBuffer[i]+nBin, temp, sizeof(double)*nBin );

			for( int i=0; i<nBin; i++ )
				tempRev[i] = temp[nBin-i-1];

			memcpy( this->domainBufferRev[i], tempRev, sizeof(double)*nBin );
			memcpy( this->domainBufferRev[i]+nBin, tempRev, sizeof(double)*nBin );

			iter++;
		}
	}

}// end copy constructor

FEL_encoder_dct&
FEL_encoder_dct::operator= ( const FEL_encoder_dct& that )
{
	if (this != &that ) // protect against invalid self-assignment
	{
		this->nBin = that.nBin;
		this->nDomain = that.nDomain;
		this->nUsedDomain = that.nUsedDomain;
		this->distanceMeasureType = that.distanceMeasureType;
		this->matchErrorThreshold = that.matchErrorThreshold;

		if( that.domainList.size() > 0 )
		{
			list<FEL_domain_core>::iterator iter = this->domainList.begin();
			this->domainList.insert( iter, that.domainList.begin(), that.domainList.end() );

			this->domainBuffer = new double*[(int)this->domainList.size()];
			this->domainBufferRev = new double*[(int)this->domainList.size()];

			iter= this->domainList.begin();
			double temp[nBin];
			double tempRev[nBin];
			for( int i=0; i<(int)this->domainList.size(); i++ )
			{
				this->domainBuffer[i] = new double[nBin*2];
				this->domainBufferRev[i] = new double[nBin*2];

				iter->getDistribution( temp );

				memcpy( this->domainBuffer[i], temp, sizeof(double)*nBin );
				memcpy( this->domainBuffer[i]+nBin, temp, sizeof(double)*nBin );

				for( int i=0; i<nBin; i++ )
					tempRev[i] = temp[nBin-i-1];

				memcpy( this->domainBufferRev[i], tempRev, sizeof(double)*nBin );
				memcpy( this->domainBufferRev[i]+nBin, tempRev, sizeof(double)*nBin );

				iter++;
			}
		}
	}
	// by convention, always return *this
	return *this;
}// end assignment operator

FEL_encoder_dct::~FEL_encoder_dct()
{
	if( domainBuffer != 0 )
	{
		for( int i=0; i<(int)this->domainList.size(); i++ )
			delete [] domainBuffer[i];
		delete [] domainBuffer;
	}
	if( domainList.size() > 0 )	domainList.clear();

}// end destructor

void
FEL_encoder_dct::init( int nbin, list<FEL_domain_core>* domainlist,
								   double threshold )
{
	nBin = nbin;
	nDomain = domainlist->size();
	nUsedDomain = 0;
	distanceMeasureType = 0;
	matchErrorThreshold = threshold;

	list<FEL_domain_core>::iterator iter = domainList.begin();
	domainList.insert( iter, domainlist->begin(), domainlist->end() );

}// end function

void
FEL_encoder_dct::organizeDomains()
{
	double domainFreqList[nBin];
	double domainGradList[nBin];
	double maxCoeffs[3];

	// Create data structure for storing reduced domains
	domainDCTList = vtkSmartPointer<vtkPoints>::New();

	#ifdef DEBUG_MODE
	FILE* domainDCTFile = fopen( "domainDCT.csv", "w" );
	#endif

	// Iterate over all domains
	vtkIdType pId = 0;
	for( list<FEL_domain_core>::iterator iter = domainList.begin();
		 iter != domainList.end();
		 iter++ )
	{
		// Get domain distribution
		iter->getDistribution( domainFreqList );

		// Compute gradient
		computeGradient( domainFreqList, domainGradList, nBin );
		#ifdef DEBUG_MODE
		fprintf( stderr, "Domain:\n" );
		for( int i=0; i<nBin-1; i++ )
			fprintf( stderr, "%g, ", domainFreqList[i] );
		fprintf( stderr, "%g\n", domainFreqList[nBin-1] );
		#endif

		// Transform using DCT and return max two coefficients
		dctTransform( domainGradList, maxCoeffs, nBin );
		#ifdef DEBUG_MODE
		fprintf( stderr, "Gradient:\n" );
		for( int i=0; i<nBin-1; i++ )
			fprintf( stderr, "%g, ", domainGradList[i] );
		fprintf( stderr, "%g\n", domainGradList[nBin-1] );
		#endif

		#ifdef DEBUG_MODE
		fprintf( domainDCTFile, "%g, %g\n", maxCoeffs[0], maxCoeffs[1] );
		#endif

		// Store coefficients in quadtree (octree)
		maxCoeffs[3] = 0;
		domainDCTList->InsertPoint( pId, maxCoeffs );
		//domainDCTList->InsertNextPoint( maxCoeffs );

		// Update range of coefficients
		if( pId == 0 )
		{
			domain2DPlaneX[0] = domain2DPlaneX[1] = maxCoeffs[0];
			domain2DPlaneY[0] = domain2DPlaneY[1] = maxCoeffs[1];
		}
		else
		{
			if( maxCoeffs[0] < domain2DPlaneX[0] ) domain2DPlaneX[0] = maxCoeffs[0];
			if( maxCoeffs[0] > domain2DPlaneX[1] ) domain2DPlaneX[1] = maxCoeffs[0];
			if( maxCoeffs[1] < domain2DPlaneY[0] ) domain2DPlaneY[0] = maxCoeffs[1];
			if( maxCoeffs[1] < domain2DPlaneY[0] ) domain2DPlaneY[1] = maxCoeffs[1];
		}

		pId ++;

	}// end for

	#ifdef DEBUG_MODE
	fclose( domainDCTFile );
	#endif

	// Compute range of DCT coefficients
	domain2DPlaneRange[0] = domain2DPlaneX[1] - domain2DPlaneX[0];
	domain2DPlaneRange[1] = domain2DPlaneY[1] - domain2DPlaneY[0];
	searchNeighborhood[0] = domain2DPlaneRange[0] * 0.01;
	searchNeighborhood[1] = domain2DPlaneRange[1] * 0.01;

	fprintf( stderr, "DCT coefficient x-range: <%g, %g>\n", domain2DPlaneX[0], domain2DPlaneX[1] );
	fprintf( stderr, "DCT coefficient y-range: <%g, %g>\n", domain2DPlaneY[0], domain2DPlaneY[1] );
	fprintf( stderr, "DCT search neighborhood: <%g, %g>\n", searchNeighborhood[0], searchNeighborhood[1] );

	// Create polydata
	fprintf( stderr, "Creating polydata ...\n" );
	domainPolyData = vtkSmartPointer<vtkPolyData>::New();
	domainPolyData->SetPoints( domainDCTList );

	//Create the octree
	fprintf( stderr, "Creating octree ...\n" );
	octree = vtkSmartPointer<vtkOctreePointLocator>::New();
	//octree = vtkSmartPointer<vtkIncrementalOctreePointLocator>::New();
	octree->SetDataSet( domainPolyData );
	octree->BuildLocator();
	fprintf( stderr, "Done\n" );
}

void
FEL_encoder_dct::computeGradient( double* fList, double* gList, int nbin )
{
	for( int i=0; i<nbin; i++ )
	{
		if( i == 0 )
			gList[i] = ( fList[i+1] - fList[i] );
		else if( i == nbin-1 )
			gList[i] = ( fList[i] - fList[i-1] );
		else
			gList[i] = ( fList[i+1] - fList[i-1] ) / 2.0;
	}
}

void
FEL_encoder_dct::dctTransform( double* gList, double* coeff, int nbin )
{
	int sz[] = {nbin, 1};
	cv::Mat data( 1, nbin, CV_64FC1 );
	cv::Mat transformedData( 1, nbin, CV_64FC1 );

	// Copy data
	for( int i=0; i<nbin; i++ )
		data.at<double>(0, i) = gList[i];
	#ifdef DEBUG_MODE
	fprintf( stderr, "Data:\n" );
	for( int i=0; i<nbin-1; i++ )
		fprintf( stderr, "%g, ", data.at<double>(0,i) );
	fprintf( stderr, "%g\n", data.at<double>(0,nBin-1) );
	#endif

	// Transform
	cv::dct( data, transformedData, 0 );
	#ifdef DEBUG_MODE
	fprintf( stderr, "DCT:\n" );
	for( int i=0; i<nbin-1; i++ )
		fprintf( stderr, "%g, ", transformedData.at<double>(0,i) );
	fprintf( stderr, "%g\n", transformedData.at<double>(0,nBin-1) );
	#endif

	// Find the maximum coefficient
	coeff[0] = transformedData.at<double>(0,0);
	for( int i=0; i<nbin; i++ )
	{
		if( transformedData.at<double>(0, i) > coeff[0] )
			coeff[0] = transformedData.at<double>(0, i);
	}

	// Find the 2nd maximum coefficient
	double diff = 100000.0;
	for( int i=0; i<nbin; i++ )
	{
		if( ( coeff[0] - transformedData.at<double>(0, i) ) < diff )
		{
			diff = ( coeff[0] - transformedData.at<double>(0, i) );
			coeff[1] = transformedData.at<double>(0, i);
		}
	}

}

void
FEL_encoder_dct::encode_Distribution_DCT( double* rangeFreqList, double* rangeLimits,
											   int searchThreshold, long nDataPoint,
					 	 	 	 	 	 	   FEL_encodedhistogram* encodedHist,
					 	 	 	 	 	 	   double errorThreshold )
{
	int optRot = 0, bestMatch = -1, nHighErrorBin = 0;
	bool optRef = false;
	double matchError = -1;
	int lowSub[3], highSub[3];
	double rangeCenter[3];
	int highErrorBinIds[nBin];
	double highErrors[nBin];
	double domainDist[nBin];
	double transformedDomain[nBin];


	#ifdef DEBUG_MODE
	double m = ITL_util<double>::Min( rangeFreqList, nBin );
	double M = ITL_util<double>::Max( rangeFreqList, nBin );
	fprintf( stderr, "in encode: Range distribution range: %g, %g\n", m, M );
	#endif

	// Get limit and center of range
	lowSub[0] = rangeLimits[0]; highSub[0] = rangeLimits[1];
	lowSub[1] = rangeLimits[2]; highSub[1] = rangeLimits[3];
	lowSub[2] = rangeLimits[4]; highSub[2] = rangeLimits[5];

	// Iterate though all the domains
	// and find the nearest one
	bestMatch = matchRangeToDomains_DCT( rangeFreqList, searchThreshold, &matchError, &optRot, &optRef );

	#ifdef DEBUG_MODE
	fprintf( stderr, "Range maps to %d-th domain with %g <%d, %d> error ...\n",
			 	     bestMatch, matchError, optRot, (int)optRef );
	#endif

	// Case 1: The range did find a match, but with poor accuracy
	// Add the high error bins to reduce error
	if( errorThreshold != 0 && matchError > errorThreshold )
	{
		//fprintf( stderr, "Storing high error bins ...\n" );

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

		// Case 1: The range did not find any match
		// Add the current range to the domain pool
		/*
		if( nHighErrorBin > 0.75*nBin )
		{
			bestMatch = correctBadEncoding( rangeFreqList, rangeLimits, &matchError, &optRot, &optRef );
			//#ifdef DEBUG_MODE
			fprintf( stderr, "After correction: Range maps to %d-th domain with %g <%d, %d> error ...\n",
					 	     	 	 	 	 	 	 	 	 bestMatch, matchError, optRot, (int)optRef );
			//#endif

			nHighErrorBin = 0;
		}
		*/
	}// end if

	// Encode range distribution
	encodedHist->setSpan( lowSub, highSub );
	encodedHist->setMatchingDomainId( bestMatch );
	encodedHist->setNumDataElement( nDataPoint );
	encodedHist->setRotation( optRot );
	encodedHist->setReflection( optRef );
	if( nHighErrorBin > 0 )
	{
		encodedHist->setNumHighErrorBin( nHighErrorBin );
		encodedHist->setHighErrorBinInfo( highErrorBinIds, highErrors );
	}
}

int
FEL_encoder_dct::matchRangeToDomains_DCT( double* rangeFreqList, int searchThreshold,
											    double *error, int* rotateAmount, bool* isReflected )
{
	int domainId = 0, bestMatchedDomainId = -1;
	double e = 0, minError = 1000000;
	double rangeDCT[3];
	double rangeGradList[nBin];
	double domainFreqList[nBin];
	double bestDomainFreqList[nBin];
	int rotAmount = 0;
	bool isRef = false;

	// Compute gradient of range
	computeGradient( rangeFreqList, rangeGradList, nBin );

	// Transform using DCT and return max two coefficients
	dctTransform( rangeGradList, rangeDCT, nBin );

	//#ifdef DEBUG_MODE
	//FILE* rangeDCTFile = fopen( "rangeDCT.csv", "a" );
	//fprintf( rangeDCTFile, "%g, %g\n", rangeDCT[0], rangeDCT[1] );
	//fclose( rangeDCTFile );
	//#endif

	// Search for domains within this range
	//vtkSmartPointer<vtkIdTypeArray> nearbyDomainIdList = vtkSmartPointer<vtkIdTypeArray>::New();
	vtkSmartPointer<vtkIdList> nearbyDomainIdList = vtkSmartPointer<vtkIdList>::New();

	//octree->FindPointsInArea( domainSearchBox, nearbyDomainIdList, true );
	octree->FindClosestNPoints( searchThreshold, rangeDCT, nearbyDomainIdList );

	//vtkIdType listSize = nearbyDomainIdList->GetNumberOfTuples();
	vtkIdType listSize = searchThreshold;
	#ifdef DEBUG_MODE
	fprintf( stderr, "Number of domains near this range: %d\n", (int)listSize );
	#endif

	// Scan through nearby domains only
	for( vtkIdType i = 0; i < listSize; i++ )
	{
		//long domainIndex = (long)nearbyDomainIdList->GetValue( i );
		vtkIdType domainInd = nearbyDomainIdList->GetId( i );
		long domainIndex = domainInd;
		//fprintf( stderr, "index: %d\n", (int)domainIndex );

		// Obtain an iterator to the domain at that position
		std::list<FEL_domain_core>::iterator domainIter = domainList.begin();
		std::advance( domainIter, domainIndex );

		// Get distribution of next domain
		domainIter->getDistribution( domainFreqList );

		// Compute match error between range and domain distributions
		//e =  FEL_util<float>::matchRangeToDomainByDistribution( distanceMeasureType,
		//														rangeFreqList, domainFreqList,
		//														nBin,
		//														&rotAmount, &isRef );

		e =  FEL_util<float>::matchRangeToDomainByDistribution_Fast( distanceMeasureType,
																	  rangeFreqList, domainIndex,
																	  domainBuffer, domainBufferRev,
																	  nBin,
																	  &rotAmount, &isRef );


		if( e < minError )
		{
			minError = e;
			(*error) = minError;
			(*rotateAmount) = rotAmount;
			(*isReflected) = isRef;
			bestMatchedDomainId = domainIndex;
			memcpy( bestDomainFreqList, domainFreqList, sizeof( double )*nBin );
		}
	}// end loop through nearby domains on the DCT space

	return bestMatchedDomainId;
}

int
FEL_encoder_dct::correctBadEncoding( double* rangeFreqList, double* rangeLimits,
										  double *error,
										  int* rotateAmount, bool* isReflected )
{
	fprintf( stderr, "Adding this range to domain list ...\n" );

	FEL_domain_core newDomain( rangeLimits, nBin, rangeFreqList );
	fprintf( stderr, "here ...\n" );
	addNewDomain( newDomain );
	fprintf( stderr, "here 1...\n" );

	// Update encoding of information of current range distribution
	(*error) = 0;
	(*rotateAmount) = 0;
	(*isReflected) = false;

	return (nDomain-1);
}

void
FEL_encoder_dct::addNewDomain( FEL_domain_core newDomain )
{
	double domainFreqList[nBin];
	double domainGradList[nBin];
	double maxCoeffs[3];

	domainList.push_back( newDomain );
	nDomain = nDomain + 1;

	fprintf( stderr, "here 3...\n" );

	// Add domain's DCT coeffs to octree
	//newDomain.getDistribution( domainFreqList );

	//computeGradient( domainFreqList, domainGradList, nBin );
	//dctTransform( domainGradList, maxCoeffs, nBin );
	//maxCoeffs[3] = 0.0;

	//fprintf( stderr, "here [%g %g %g]...\n", maxCoeffs[0], maxCoeffs[1], maxCoeffs[2] );
	//int id = octree->IsInsertedPoint( maxCoeffs );
	//if( id == -1 )
	//	octree->InsertNextPoint( maxCoeffs );
	//else
	//	fprintf( stderr, "Point already in octree %d ...\n", id );

	//#ifdef DEBUG_MODE
	fprintf( stderr, "Total domain count: %d\n", nDomain );
	//#endif
}

int
FEL_encoder_dct::getDomainCount()
{
	return nDomain;
}
