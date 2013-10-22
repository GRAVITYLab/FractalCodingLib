/*
 * FEL_encoder_incremental.h
 *
 *  Created on: Sep 10, 2012
 *      Author: abon
 */

#ifndef FEL_ENCODER_DCT_H_
#define FEL_ENCODER_DCT_H_

#include <list>

#include <cv.h>
#include <highgui.h>

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
#include <vtkIncrementalOctreePointLocator.h>
#include <vtkLookupTable.h>
#include <vtkIdTypeArray.h>
#include <vtkIdList.h>
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

class FEL_encoder_dct
{
	int nBin, nDomain, nUsedDomain;
	int distanceMeasureType;
	double matchErrorThreshold;

	list<FEL_domain_core> domainList;
	double** domainBuffer;
	double** domainBufferRev;

	vtkSmartPointer<vtkPoints> domainDCTList;
	vtkSmartPointer<vtkPolyData> domainPolyData;
	vtkSmartPointer<vtkOctreePointLocator> octree;
	//vtkSmartPointer<vtkIncrementalOctreePointLocator> octree;

	double domain2DPlaneX[2], domain2DPlaneY[2], domain2DPlaneRange[2];
	double searchNeighborhood[2];

public:

	// Constructor
	FEL_encoder_dct();
	FEL_encoder_dct( int nbin, list<FEL_domain_core>* domainlist, double threshold );

	FEL_encoder_dct& operator= ( const FEL_encoder_dct& that );
	FEL_encoder_dct( const FEL_encoder_dct& that );

	~FEL_encoder_dct();

	void init( int nbin, list<FEL_domain_core>* domainlist, double threshold );

	void organizeDomains();
	void computeGradient( double* fList, double* gList, int nbin );
	void dctTransform( double* gList, double* coeffs, int nbin );


	void encode_Distribution_DCT( double* rangeFreqList, double* rangeLimits,
									  int searchThreshold, long nDataPoint,
									  FEL_encodedhistogram* encodedHist,
									  double errorThreshold );

	int matchRangeToDomains_DCT( double* rangeFreqList, int searchThreshold,
		    						double *error, int* rotateAmount, bool* isReflected );

	void addNewDomain( FEL_domain_core newDomain );
	int correctBadEncoding( double* rangeFreqList, double* rangeLimits,
			  	  	  	  	   double *error,
			  	  	  	  	   int* rotateAmount, bool* isReflected );

	int getDomainCount();

};

#endif /* FEL_ENCODER_DCT_H_ */
