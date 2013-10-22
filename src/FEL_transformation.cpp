/*
 * FEL_transformation.cpp
 *
 *  Created on: Feb 15, 2012
 *      Author: abon
 */

#include "FEL_transformation.h"

FEL_transformation::FEL_transformation()
{
}

FEL_transformation::FEL_transformation( int nrot )
{
	nRotation = nrot;
	nReflection = 2;

	rotMatrixArray = new MATRIX3[nRotation];
	refMatrixArray = new MATRIX3[nReflection];

	initRotationMatrices();
	initReflectionMatrices();

}// End constructor

void
FEL_transformation::initTransformationMatrices( int nrot )
{
	nRotation = nrot;
	nReflection = 2;

	rotMatrixArray = new MATRIX3[nRotation];
	refMatrixArray = new MATRIX3[nReflection];

	initRotationMatrices();
	initReflectionMatrices();

}

void
FEL_transformation::initRotationMatrices()
{
	float angleStep = 360.0f / nRotation;
	for( int i=0; i<nRotation; i++ )
	{
		float theta = (pi/180.0f) * (i * angleStep);
		VECTOR3 r0( cos(theta), -sin(theta), 0 );
		VECTOR3 r1( sin(theta), cos(theta), 0 );
		VECTOR3 r2( 0, 0, 1 );

		rotMatrixArray[i][0] = r0;
		rotMatrixArray[i][1] = r1;
		rotMatrixArray[i][2] = r2;

	}
}// End function

void
FEL_transformation::initReflectionMatrices()
{
	VECTOR3 r0( 1, 0, 0 );
	VECTOR3 r1( 0, -1, 0 );
	VECTOR3 r2( 0, 0, 1 );

	refMatrixArray[0][0] = r0;
	refMatrixArray[0][1] = r1;
	refMatrixArray[0][2] = r2;

	r0.Set( -1, 0, 0 );
	r1.Set( 0, 1, 0 );
	r2.Set( 0, 0, 1 );

	refMatrixArray[1][0] = r0;
	refMatrixArray[1][1] = r1;
	refMatrixArray[1][2] = r2;

}// End function

void
FEL_transformation::applyTransformation( int type, int id, VECTOR3* data, int npoint )
{
	for( int iP=0; iP<npoint; iP++ )
	{
		VECTOR3 tmp( data[iP] );

		if( type == ROTATION )
			tmp = FEL_transformation::mult( rotMatrixArray[id], tmp );
		else if( type == REFLECTION )
			tmp = FEL_transformation::mult( refMatrixArray[id], tmp );

		data[iP].Set( tmp(0), tmp(1), tmp(2) );
	}
}

void
FEL_transformation::applyTransformationToDistribution( int type, int nStep, float *freqList, int nBin )
{
	if( type == ROTATION )
	{
		//float tmp1[nBin];
		//float tmp2[nBin];
		//int part = nBin - nStep;
		//memcpy( tmp1, freqList, sizeof(float)*part );
		//memcpy( tmp2, freqList+part, sizeof(float)*nStep );

		//memcpy( freqList, tmp2, sizeof(float)*nStep );
		//memcpy( freqList+nStep, tmp1, sizeof(float)*part );
	}

}

// return m0 * v0
VECTOR3
FEL_transformation::mult( MATRIX3 m0, VECTOR3 v0 )
{
	VECTOR3 result;

	result[0] = dot(m0(0),v0);
	result[1] = dot(m0(1),v0);
	result[2] = dot(m0(2),v0);

	return(result);
}

FEL_transformation::~FEL_transformation()
{
	delete [] rotMatrixArray;
	delete [] refMatrixArray;
}// End destructor

