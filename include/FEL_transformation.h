/*
 * FEL_transformation.h
 *
 *  Created on: Feb 15, 2012
 *      Author: abon
 */

#ifndef FEL_TRANSFORMATION_H_
#define FEL_TRANSFORMATION_H_

#include "ITL_header.h"
#include "ITL_vectormatrix.h"

class FEL_transformation
{
	enum transformTypes { ROTATION = 0,
						  REFLECTION = 1 };

	MATRIX3 *rotMatrixArray;
	MATRIX3 *refMatrixArray;

public:

	int nRotation;
	int nReflection;

	FEL_transformation();
	FEL_transformation( int nrot );

	void initTransformationMatrices( int nrot );
	void initRotationMatrices();
	void initReflectionMatrices();

	void applyTransformation( int type, int id, VECTOR3* data, int npoint );

	void applyTransformationToDistribution( int type, int nStep, float *freqList, int nBin );

	static VECTOR3 mult( MATRIX3 m0, VECTOR3 v0 );

	~FEL_transformation();
};

#endif
/* FEL_TRANSFORMATION_H_ */
