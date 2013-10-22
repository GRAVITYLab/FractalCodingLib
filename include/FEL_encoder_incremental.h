/*
 * FEL_encoder_incremental.h
 *
 *  Created on: Sep 10, 2012
 *      Author: abon
 */

#ifndef FEL_ENCODER_INCREMENTAL_H_
#define FEL_ENCODER_INCREMENTAL_H_

#include <list>

//#include "ITL_distcomputer.h"
//#include "ITL_distribution.h"
#include "ITL_histogram.h"
#include "ITL_entropycore.h"
#include "ITL_field_regular.h"

#include "FEL_util.h"
#include "FEL_segmenttree2.h"
#include "FEL_domain_core.h"
#include "FEL_encodedhistogram.h"

class FEL_encoder_incremental
{
	enum transformTypes { ROTATION = 0,
						  REFLECTION = 1 };

	int nBin, nDomain, nUsedDomain;
	float minDomainEntropy, maxDomainEntropy;
	int distanceMeasureType;
	double matchErrorThreshold;

	int domainUtilizationList[50000];
	list<FEL_domain_core> domainList;
	FEL_segmenttree2 domainTree;

public:

	// Constructor
	FEL_encoder_incremental();
	FEL_encoder_incremental( int nbin, list<FEL_domain_core>* domainlist, double threshold );

	FEL_encoder_incremental& operator= ( const FEL_encoder_incremental& that );
	FEL_encoder_incremental( const FEL_encoder_incremental& that );

	~FEL_encoder_incremental();

	void init( int nbin, list<FEL_domain_core>* domainlist, double threshold );

	void computeDomainEntropies();
	void sortDomains_Entropybased();
	//void pruneDomains_Entropybased();
	//void pruneDomains_Entropybased2();
	void organizeDomains(  bool isDomainPruningOn );
	FEL_segmenttree2* getDomainTree();

	int encode_Distribution( double* rangeFreqList, long nDataPoint,
							  FEL_encodedhistogram* encodedHist,
							  double errorThreshold, bool allZeroFlag = false );

	int encode_Distribution_2( double* spanFreqList,
								  double* residualFreqList,
								  long nDataPoint,
								  FEL_encodedhistogram* encodedHist );

	int encode_Distribution_2( double* spanFreqList,
								  double* residualFreqList,
								  long nDataPoint,
		 	 	 	 	  	  	  int* optTemplateId,
		 	 	 	 	  	  	  int* optRotation,
		 	 	 	 	  	  	  bool* reflectionFlag );

	int findBestMatch_peakbased( double* freqList,
									 double* errorList,
									 int* rotateAmount,
									 bool* isReflected );

	int matchRangeToDomains2( double* rangeFreqList,
								 double *error,
								 int* rotateAmount, bool* isReflected );
	int matchRangeToDomains_Entropybased( double* rangeFreqList,
											   double *error,
											   int* rotateAmount, bool* isReflected );

	void transformDistribution( double* freqList, int nRot, bool isRef );

	void updateDomainSet();

	void addNewDomain( FEL_domain_core newDomain );
	void correctBadEncoding(  FEL_encodedhistogram* encodedHist,
		 						 double* rangeFreqList, int nData );

	int getDomainCount();
	int getActiveDomainCount();

	void incrementDomainUtilization( int id );
	void decrementDomainUtilization( int id );

};

#endif /* FEL_ENCODER_INCREMENTAL_H_ */
