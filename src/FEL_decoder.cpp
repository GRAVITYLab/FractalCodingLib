/*
 * FEL_decoder.cpp
 *
 *  Created on: Mar 20, 2013
 *      Author: abon
 */

#include "FEL_decoder.h"

FEL_decoder::FEL_decoder()
{

}

FEL_decoder::FEL_decoder( int nbin, list<FEL_domain_core>* domainlist )
{
	nBin = nbin;
	nDomain = domainlist->size();

	list<FEL_domain_core>::iterator iter = domainList.begin();
	domainList.insert( iter, domainlist->begin(), domainlist->end() );
}

FEL_decoder::FEL_decoder( const FEL_decoder& that )
{
	this->nBin = that.nBin;
	this->nDomain = that.nDomain;

	if( that.domainList.size() > 0 )
	{
		list<FEL_domain_core>::iterator iter = this->domainList.begin();
		this->domainList.insert( iter, that.domainList.begin(), that.domainList.end() );
	}
}

FEL_decoder&
FEL_decoder::operator= ( const FEL_decoder& that )
{
	if (this != &that ) // protect against invalid self-assignment
	{
		this->nBin = that.nBin;
		this->nDomain = that.nDomain;

		if( that.domainList.size() > 0 )
		{
			list<FEL_domain_core>::iterator iter = this->domainList.begin();
			this->domainList.insert( iter, that.domainList.begin(), that.domainList.end() );
		}
	}
	// by convention, always return *this
	return *this;
}

FEL_decoder::~FEL_decoder()
{

}

void
FEL_decoder::init( int nbin, list<FEL_domain_core>* domainlist )
{
	nBin = nbin;
	nDomain = domainlist->size();

	list<FEL_domain_core>::iterator iter = domainList.begin();
	domainList.insert( iter, domainlist->begin(), domainlist->end() );

}// end function

long
FEL_decoder::decode_Distribution( double* decodedFreqList,
									  FEL_encodedhistogram* encodedHist,
									  bool normFlag )
{
	double freqList[nBin];

	// Get information from encoder
	int matchingDomainId = encodedHist->getMatchingDomainId();
	long nDataElement = encodedHist->getNumDataElement();
	int nHighErrorBin = encodedHist->getNumHighErrorBin();
	int optRot = encodedHist->getRotation();
	bool optRef = encodedHist->getReflection();

	#ifdef DEBUG_MODE
	fprintf( stderr, "In decoder: matching domain id: %d\n", matchingDomainId );
	fprintf( stderr, "In decoder: number of elements: %d\n", nDataElement );
	fprintf( stderr, "In decoder: rotation and reflection: %d %d\n", optRot, optRef );
	fprintf( stderr, "In decoder: Number of high error bins: %d\n", nHighErrorBin );
	fprintf( stderr, "In decoder: matching domain distribution:\n" );
	#endif

	// Get domain distribution with specified ID
	list<FEL_domain_core>::iterator domainIter = domainList.begin();
	for( int iD=0; iD<matchingDomainId; iD++ )
		domainIter++;
	domainIter->getDistribution( freqList );

	#ifdef DEBUG_MODE
	for( int i = 0; i<nBin-1; i++ )
		fprintf( stderr, "%g, ", freqList[i] );
	fprintf( stderr, "%g\n", freqList[nBin-1] );
	#endif

	// Transform domain to match range
	FEL_util<double>::transformDistribution( freqList, nBin, optRot, optRef );
	#ifdef DEBUG_MODE
	fprintf( stderr, "In decoder: transformed domain distribution:\n" );
	for( int i = 0; i<nBin-1; i++ )
		fprintf( stderr, "%g, ", freqList[i] );
	fprintf( stderr, "%g\n", freqList[nBin-1] );
	#endif

	// Further refine the high error bins
	if( nHighErrorBin > 0 )
	{
		int highErrorBinIds[nHighErrorBin];
		double highErrors[nHighErrorBin];

		//fprintf( stderr, "In decoder: here to correct high error\n" );

		// Get stored error information
		encodedHist->getHighErrorBinInfo( highErrorBinIds, highErrors );

		#ifdef DEBUG_MODE
		fprintf( stderr, "In decoder: Error Info:\n" );
		for( int i = 0; i<nHighErrorBin; i++ )
			fprintf( stderr, "%d: %g\n", highErrorBinIds[i], highErrors[i] );
		#endif

		// Use information to improve decoding
		correctHighErrorBins( freqList, nHighErrorBin, highErrorBinIds, highErrors );
		#ifdef DEBUG_MODE
		fprintf( stderr, "In decoder: error corrected decoded distribution:\n" );
		for( int i = 0; i<nBin-1; i++ )
			fprintf( stderr, "%g, ", freqList[i] );
		fprintf( stderr, "%g\n", freqList[nBin-1] );
		#endif

		// Renormalize
		double sumFreq = ITL_util<double>::sum( freqList, nBin );
		ITL_util<double>::divideArrayScalar( freqList, sumFreq, nBin );
		#ifdef DEBUG_MODE
		fprintf( stderr, "In decoder: re-normalized decoded distribution:\n" );
		for( int i = 0; i<nBin-1; i++ )
			fprintf( stderr, "%g, ", freqList[i] );
		fprintf( stderr, "%g\n", freqList[nBin-1] );
		#endif

	}

	//fprintf( stderr, "In decoder: here 1\n" );
	// Multiply number of elements with normalized distribution
	if( nDataElement != 0 && normFlag == false )
	{
		for( int i=0; i<nBin; i++ )
			freqList[i] = freqList[i] * nDataElement;
	}

	// Copy decoded distribution
	//fprintf( stderr, "In decoder: here 2\n" );
	memcpy( decodedFreqList, freqList, sizeof(double)*nBin );

	return nDataElement;

}// end function

void
FEL_decoder::decode_Distribution_2( double* decodedFreqList,
										 double* residueHist,
									  	 FEL_encodedhistogram* encodedHist,
									  	 bool normFlag )
{
	double templateHist[nBin];
	double shiftedtemplateHist[nBin];

	// Get information from encoder
	int matchingDomainId = encodedHist->getMatchingDomainId();
	long nDataElement = encodedHist->getNumDataElement();
	//int nHighErrorBin = encodedHist->getNumHighErrorBin();
	int optRot = encodedHist->getRotation();
	//bool optRef = encodedHist->getReflection();
	//fprintf( stderr, "template id: %d\n", matchingDomainId );
	//fprintf( stderr, "shift: %d\n", optRot );

	// Get domain distribution with specified ID
	list<FEL_domain_core>::iterator domainIter = domainList.begin();
	for( int iD=0; iD<matchingDomainId; iD++ )
		domainIter++;
	domainIter->getDistribution( templateHist );
	//fprintf( stderr, "here\n" );

	// Transform domain to match range
	//FEL_util<double>::shiftDistribution_zerofill( templateHist, shiftedtemplateHist, nBin, optRot );
	FEL_util<double>::shiftDistribution_rotate( templateHist, shiftedtemplateHist, nBin, optRot );
	//fprintf( stderr, "here 2\n" );

	// Adjust errors
	for( int i=0; i<nBin; i++ )
		decodedFreqList[i] = shiftedtemplateHist[i] + residueHist[i];
	//fprintf( stderr, "here 3\n" );

	//fprintf( stderr, "In decoder: here 1\n" );
	// Multiply number of elements with normalized distribution
	if( nDataElement != 0 && normFlag == false )
	{
		for( int i=0; i<nBin; i++ )
			decodedFreqList[i] = decodedFreqList[i] * nDataElement;
	}

}// end function

void
FEL_decoder::decode_Distribution_2( double* decodedFreqList,
										int matchingDomainId,
										int optRot,
										bool optRef,
										int nDataElement,
										double* residueHist,
									  	bool normFlag )
{
	double templateHist[nBin];
	double shiftedtemplateHist[nBin];

	// Get domain distribution with specified ID
	list<FEL_domain_core>::iterator domainIter = domainList.begin();
	for( int iD=0; iD<matchingDomainId; iD++ )
		domainIter++;
	domainIter->getDistribution( templateHist );
	//fprintf( stderr, "here\n" );

	// Transform domain to match range
	//FEL_util<double>::shiftDistribution_zerofill( templateHist, shiftedtemplateHist, nBin, optRot );
	FEL_util<double>::shiftDistribution_rotate( templateHist, shiftedtemplateHist, nBin, optRot );
	//fprintf( stderr, "here 2\n" );

	// Adjust errors
	for( int i=0; i<nBin; i++ )
		decodedFreqList[i] = shiftedtemplateHist[i] + residueHist[i];
	//fprintf( stderr, "here 3\n" );

	//fprintf( stderr, "In decoder: here 1\n" );
	// Multiply number of elements with normalized distribution
	if( nDataElement != 0 && normFlag == false )
	{
		for( int i=0; i<nBin; i++ )
			decodedFreqList[i] = decodedFreqList[i] * nDataElement;
	}

}// end function





void
FEL_decoder::correctHighErrorBins( double* freqList, int nHighErrorBin,
											   int* highErrorBinId, double* highError )
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
