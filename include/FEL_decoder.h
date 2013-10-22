/*
 * FEL_decoder.h
 *
 *  Created on: Mar 20, 2013
 *      Author: abon
 */

#ifndef FEL_DECODER_H_
#define FEL_DECODER_H_

#include "FEL_util.h"
#include "FEL_domain_core.h"
#include "FEL_encodedhistogram.h"

class FEL_decoder
{
private:

	int nBin, nDomain;
	list<FEL_domain_core> domainList;

public:

	// Constructor
	FEL_decoder();
	FEL_decoder( int nbin, list<FEL_domain_core>* domainlist );

	FEL_decoder& operator= ( const FEL_decoder& that );
	FEL_decoder( const FEL_decoder& that );

	~FEL_decoder();

	void init( int nbin, list<FEL_domain_core>* domainlist );

	void transformDistribution( double* freqList, int nRot, bool isRef );

	long decode_Distribution( double* decodedFreqList, FEL_encodedhistogram* encodedHist, bool normFlag = false );
	void decode_Distribution_2( double* decodedFreqList, double* residueHist,
								    FEL_encodedhistogram* encodedHist, bool normFlag );
	void decode_Distribution_2( double* decodedFreqList,
									int matchingDomainId, int optRot, bool optRef,
									int nDataElement, double* residueHist, bool normFlag );


	void correctHighErrorBins( double* freqList, int nHighErrorBin, int* highErrorBinId, double* highError );
};

#endif
/* FEL_DECODER_H_ */
