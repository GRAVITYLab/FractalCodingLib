/*
 * FEL_encoder_localsearch.h
 *
 *  Created on: Mar 20, 2013
 *      Author: abon
 */

#ifndef FEL_ENCODER_LOCALSEARCH_H_
#define FEL_ENCODER_LOCALSEARCH_H_

#include <list>

#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkProperty.h>
#include <vtkDataSetMapper.h>
#include <vtkPolyDataMapper.h>
#include <vtkCubeSource.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkPointSource.h>
#include <vtkPointData.h>
#include <vtkOctreePointLocator.h>
#include <vtkLookupTable.h>
#include <vtkIdTypeArray.h>
#include <vtkCommand.h>
#include <vtkDenseArray.h>
#include "vtkSmartPointer.h"
#include "vtkByteSwap.h"
#include "vtkDoubleArray.h"

#include "ITL_distcomputer.h"
#include "ITL_distribution.h"
#include "ITL_histogram.h"
#include "ITL_entropycore.h"
#include "ITL_field_regular.h"

#include "FEL_util.h"
#include "FEL_domain_core.h"
#include "FEL_encodedhistogram.h"

class FEL_encoder_localsearch
{
	int nBin, nDomain, nUsedDomain;
	float minDomainEntropy, maxDomainEntropy;
	int distanceMeasureType;
	double matchErrorThreshold;

	list<FEL_domain_core> domainList;

	vtkSmartPointer<vtkPoints> domainCenterList;
	vtkSmartPointer<vtkPolyData> domainPolyData;
	vtkSmartPointer<vtkOctreePointLocator> octree;

public:

	// Constructor
	FEL_encoder_localsearch();
	FEL_encoder_localsearch( int nbin, list<FEL_domain_core>* domainlist, double threshold );

	FEL_encoder_localsearch& operator= ( const FEL_encoder_localsearch& that );
	FEL_encoder_localsearch( const FEL_encoder_localsearch& that );

	~FEL_encoder_localsearch();

	void init( int nbin, list<FEL_domain_core>* domainlist, double threshold );

	void organizeDomains();

	void encode_Distribution_LocalSearch( double* rangeFreqList, double* rangeCenter,
											  double distanceThreshold, long nDataPoint,
							  	  	  	  	  FEL_encodedhistogram* encodedHist,
							  	  	  	  	  double errorThreshold, bool allZeroFlag = false );

	int matchRangeToDomains_LocalSearch( double* rangeFreqList, double* rangeCenter,
			  	  	  	  	  	  	  	  	  double distanceThreshold, double *error,
											  int* rotateAmount, bool* isReflected );

	void addNewDomain( FEL_domain_core newDomain );
	int correctBadEncoding(  double* rangeFreqList, double* rangeLimits,
								double *error, 	int* rotateAmount, bool* isReflected );

	int getDomainCount();
};

#endif
/* FEL_ENCODER_LOCALSEARCH_H_ */
