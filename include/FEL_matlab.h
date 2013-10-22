/*
 * FEL_matlab.h
 *
 *  Created on: Mar 5, 2012
 *      Author: abon
 */

#ifndef FEL_MATLAB_H_
#define FEL_MATLAB_H_

#include <iostream>
#include <cstring>
#include <cstdio>
#include "engine.h"

#define NMAX 25000

class FEL_matlab
{
	Engine *engine;

	mxArray *mxNBin, *mxNTotalDomain, *mxNTotalRange,
			*mxNUnmatchedRange, *mxUnmatchedRangeFlagList,
			*mxErrorList, *mxMuErrorList, *mxSdErrorList,
			*mxMatchedDomainList,
			*mxDomainMatchFreqList,
			*mxRangeLimitList,
			*mxDomEntropyList, *mxRangeEntropyList;

public:

	int nBin, nDomain, nRange, nUnmatchedRange;
	double errorList[NMAX];
	double muErrorList[NMAX];
	double sdErrorList[NMAX];
	double matchedDomainIdList[NMAX];
	double rangeSizeList[NMAX];
	double domainMatchFreqList[NMAX];
	double domainEntropyList[NMAX];
	double rangeEntropyList[NMAX];
	double rangeLimitList[7*NMAX];
	double unMatchedRangeFlagList[NMAX];

	FEL_matlab();
	~FEL_matlab();

	void initEngine();
	void closeEngine();

	void createMatlabVariables();

	void displayDomainUtilization_Matlab();
	void displayMatchingPairEntropyLevel_Matlab();
	void displayMatchErrorDistribution_Matlab();
	void displayBlockwiseMatchErrorPlot_Matlab();
	void displayUnMatchedRanges_Matlab();

	void displayDomainEntropies_Matlab();


};

#endif
/* FEL_MATLAB_H_ */
