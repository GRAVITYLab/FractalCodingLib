/*
 * FEL_encodedhistogram.cpp
 *
 *  Created on: Sep 9, 2012
 *      Author: abon
 */

#include "FEL_encodedhistogram.h"

FEL_encodedhistogram::FEL_encodedhistogram()
{
	// Zero-set range
	memset( low, 0, sizeof(int)*3 );
	memset( high, 0, sizeof(int)*3 );

	spanId = 0;

	matchingDomainId = 0;
	nElement = 0;
	nHighErrorBin = 0;
	nRotation = 0;
	reflectionFlag = false;

	highErrorBinIdList = NULL;
	highErrorList = NULL;
}

FEL_encodedhistogram::FEL_encodedhistogram( int domainId, int nElem, int nHighErrBin,
									int* highErrBinId, double* highErrList )
{
	// Zero-set range
	memset( low, 0, sizeof(int)*3 );
	memset( high, 0, sizeof(int)*3 );

	spanId = 0;

	matchingDomainId = domainId;
	nElement = nElem;
	nHighErrorBin = nHighErrBin;
	nRotation = 0;
	reflectionFlag = false;

	highErrorBinIdList = NULL;
	highErrorList = NULL;

	if( nHighErrorBin > 0 )
	{
		highErrorBinIdList = new int[nHighErrorBin];
		highErrorList = new double[nHighErrorBin];

		memcpy( highErrorBinIdList, highErrBinId, sizeof(int)*nHighErrorBin );
		memcpy( highErrorList, highErrList, sizeof(double)*nHighErrorBin );
	}

}// End constructor

FEL_encodedhistogram::FEL_encodedhistogram( int* l, int* h,
		 										  int domainId,
		 										  int nElem, int nHighErrBin,
		 										  int* highErrBinId, double* highErrList )
{
	// Set range
	memcpy( low, l, sizeof(int)*3 );
	memcpy( high, l, sizeof(int)*3 );

	spanId = 0;

	matchingDomainId = domainId;
	nElement = nElem;
	nHighErrorBin = nHighErrBin;
	nRotation = 0;
	reflectionFlag = false;

	highErrorBinIdList = NULL;
	highErrorList = NULL;

	if( nHighErrorBin > 0 )
	{
		highErrorBinIdList = new int[nHighErrorBin];
		highErrorList = new double[nHighErrorBin];

		memcpy( highErrorBinIdList, highErrBinId, sizeof(int)*nHighErrorBin );
		memcpy( highErrorList, highErrList, sizeof(double)*nHighErrorBin );
	}

}// End constructor


FEL_encodedhistogram::FEL_encodedhistogram( const FEL_encodedhistogram& that )
{
	// Set range
	memcpy( this->low, that.low, sizeof(int)*3 );
	memcpy( this->high, that.high, sizeof(int)*3 );

	this->spanId = that.spanId;

	this->matchingDomainId = that.matchingDomainId;
	this->nElement = that.nElement;
	this->nHighErrorBin = that.nHighErrorBin;
	this->nRotation = that.nRotation;
	this->reflectionFlag = that.reflectionFlag;

	this->highErrorBinIdList = NULL;
	this->highErrorList = NULL;

	if( this->nHighErrorBin > 0 )
	{
		this->highErrorBinIdList = new int[nHighErrorBin];
		this->highErrorList = new double[nHighErrorBin];

		memcpy( this->highErrorBinIdList, that.highErrorBinIdList, sizeof(int)*nHighErrorBin );
		memcpy( this->highErrorList, that.highErrorList, sizeof(double)*nHighErrorBin );
	}

}

FEL_encodedhistogram&
FEL_encodedhistogram::operator= ( const FEL_encodedhistogram& that )
{
	if (this != &that ) // protect against invalid self-assignment
	{
		// Set range
		memcpy( this->low, that.low, sizeof(int)*3 );
		memcpy( this->high, that.high, sizeof(int)*3 );

		this->spanId = that.spanId;

		this->matchingDomainId = that.matchingDomainId;
		this->nElement = that.nElement;
		this->nHighErrorBin = that.nHighErrorBin;
		this->nRotation = that.nRotation;
		this->reflectionFlag = that.reflectionFlag;

		this->highErrorBinIdList = NULL;
		this->highErrorList = NULL;

		if( this->nHighErrorBin > 0 )
		{
			this->highErrorBinIdList = new int[nHighErrorBin];
			this->highErrorList = new double[nHighErrorBin];

			memcpy( this->highErrorBinIdList, that.highErrorBinIdList, sizeof(int)*nHighErrorBin );
			memcpy( this->highErrorList, that.highErrorList, sizeof(double)*nHighErrorBin );
		}
	}

	// by convention, always return *this
	return *this;
}

FEL_encodedhistogram::~FEL_encodedhistogram()
{
	if( highErrorBinIdList != NULL )	delete [] highErrorBinIdList;
	if( highErrorList != NULL )			delete [] highErrorList;
}

void
FEL_encodedhistogram::getSpan( int* l, int* h )
{
	//memcpy( l, this->low, sizeof(int)*3 );
	//memcpy( h, this->high, sizeof(int)*3 );
	for( int i=0; i<3; i++ )
	{
		l[i] = low[i];
		h[i] = high[i];
	}
}

int
FEL_encodedhistogram::getSpan2( int* l, int* h )
{
	memcpy( l, this->low, sizeof(int)*3 );
	memcpy( h, this->high, sizeof(int)*3 );
	return spanId;
}

int
FEL_encodedhistogram::getSpanId()
{
	return spanId;
}

int
FEL_encodedhistogram::getMatchingDomainId()
{
	return matchingDomainId;
}

long
FEL_encodedhistogram::getNumDataElement()
{
	return nElement;
}

double
FEL_encodedhistogram::getEncodingError()
{
	return encodingError;
}

int
FEL_encodedhistogram::getNumHighErrorBin()
{
	return nHighErrorBin;
}

int
FEL_encodedhistogram::getRotation()
{
	return nRotation;
}

bool
FEL_encodedhistogram::getReflection()
{
	return reflectionFlag;
}

void
FEL_encodedhistogram::getHighErrorBinInfo( int* binIds, double* highErrors )
{
	assert( binIds != NULL );
	assert( highErrors != NULL );
	assert( highErrorBinIdList != NULL );
	assert( highErrorList != NULL );

	if( nHighErrorBin == 0 )
		return;

	//fprintf( stderr, "nHigh: %d\n", nHighErrorBin );

	memcpy( binIds, highErrorBinIdList, sizeof(int)*nHighErrorBin );
	memcpy( highErrors, highErrorList, sizeof(double)*nHighErrorBin );

	//fprintf( stderr, "nHigh: %d\n", nHighErrorBin );
}

void
FEL_encodedhistogram::getHighErrorBinIds( int* binIds )
{
	assert( binIds != NULL );
	assert( highErrorBinIdList != NULL );

	if( nHighErrorBin == 0 )
		return;

	memcpy( binIds, highErrorBinIdList, sizeof(int)*nHighErrorBin );
}

void
FEL_encodedhistogram::getHighErrorBinFreqs( double* highErrors )
{
	assert( highErrors != NULL );
	assert( highErrorList != NULL );

	if( nHighErrorBin == 0 )
		return;

	memcpy( highErrors, highErrorList, sizeof(double)*nHighErrorBin );
}

void
FEL_encodedhistogram::setSpan( int* l, int* h )
{
	memcpy( low, l, sizeof(int)*3 );
	memcpy( high, h, sizeof(int)*3 );
}

void
FEL_encodedhistogram::setSpan( int* l, int* h, int id )
{
	//memcpy( low, l, sizeof(int)*3 );
	//memcpy( high, h, sizeof(int)*3 );
	for( int i=0; i<3; i++ )
	{
		low[i] = l[i];
		high[i] = h[i];
	}
	spanId = id;
}

void
FEL_encodedhistogram::setSpanId( int id )
{
	spanId = id;
}

void
FEL_encodedhistogram::setMatchingDomainId( int n )
{
	matchingDomainId = n;
}

void
FEL_encodedhistogram::setNumDataElement( long n )
{
	nElement = n;
}

void
FEL_encodedhistogram::setEncodingError( double d )
{
	encodingError = d;
}


void
FEL_encodedhistogram::setNumHighErrorBin( int n )
{
	nHighErrorBin = n;
}

void
FEL_encodedhistogram::setRotation( int nrot )
{
	nRotation = nrot;
}

void
FEL_encodedhistogram::setReflection( bool isref )
{
	reflectionFlag = isref;
}

void
FEL_encodedhistogram::setHighErrorBinInfo( int* binIds, double* highErrors )
{
	if( nHighErrorBin == 0 )
		return;

	if( highErrorBinIdList == NULL ) highErrorBinIdList = new int[nHighErrorBin];
	if( highErrorList == NULL ) highErrorList = new double[nHighErrorBin];

	memcpy( highErrorBinIdList, binIds, sizeof(int)*nHighErrorBin );
	memcpy( highErrorList, highErrors, sizeof(double)*nHighErrorBin );
}

void
FEL_encodedhistogram::setHighErrorBinIds( int* binIds )
{
	if( nHighErrorBin == 0 )
		return;

	if( highErrorBinIdList == NULL ) highErrorBinIdList = new int[nHighErrorBin];

	memcpy( highErrorBinIdList, binIds, sizeof(int)*nHighErrorBin );
}

void
FEL_encodedhistogram::setHighErrorBinFreqs( double* highErrors )
{
	if( nHighErrorBin == 0 )
		return;
	if( highErrorList == NULL ) highErrorList = new double[nHighErrorBin];

	memcpy( highErrorList, highErrors, sizeof(double)*nHighErrorBin );
}

void
FEL_encodedhistogram::read( FILE* pFile )
{
	fread( low, sizeof(int), 3, pFile );
	fread( high, sizeof(int), 3, pFile );

	fread( &matchingDomainId, sizeof(int), 1, pFile );
	fread( &nRotation, sizeof(int), 1, pFile );
	fread( &reflectionFlag, sizeof(bool), 1, pFile );
	fread( &nElement, sizeof(long), 1, pFile );
	//fread( &nHighErrorBin, sizeof(int), 1, pFile );
	//if( nHighErrorBin > 0 )
	//{
	//	highErrorBinIdList = new int[nHighErrorBin];
	//	highErrorList = new double[nHighErrorBin];

	//	fread( highErrorBinIdList, sizeof(int), nHighErrorBin, pFile );
	//	fread( highErrorList, sizeof(double), nHighErrorBin, pFile );
	//}
}


void
FEL_encodedhistogram::write( FILE* pFile )
{
	fwrite( low, sizeof(int), 3, pFile );
	fwrite( high, sizeof(int), 3, pFile );

	fwrite( &matchingDomainId, sizeof(int), 1, pFile );
	fwrite( &nRotation, sizeof(int), 1, pFile );
	fwrite( &reflectionFlag, sizeof(bool), 1, pFile );
	fwrite( &nElement, sizeof(long), 1, pFile );
	//fwrite( &nHighErrorBin, sizeof(int), 1, pFile );
	//if( nHighErrorBin > 0 )
	//{
	//	fwrite( highErrorBinIdList, sizeof(int), nHighErrorBin, pFile );
	//	fwrite( highErrorList, sizeof(double), nHighErrorBin, pFile );
	//}
}

void
FEL_encodedhistogram::read2( FILE* pFile )
{
	fread( &spanId, sizeof(int), 1, pFile );

	fread( &matchingDomainId, sizeof(int), 1, pFile );
	fread( &nRotation, sizeof(int), 1, pFile );
	fread( &reflectionFlag, sizeof(bool), 1, pFile );
	fread( &nHighErrorBin, sizeof(int), 1, pFile );
	if( nHighErrorBin > 0 )
	{
		highErrorBinIdList = new int[nHighErrorBin];
		highErrorList = new double[nHighErrorBin];

		fread( highErrorBinIdList, sizeof(int), nHighErrorBin, pFile );
		fread( highErrorList, sizeof(double), nHighErrorBin, pFile );
	}
}

void
FEL_encodedhistogram::write2( FILE* pFile  )
{
	fwrite( &spanId, sizeof(int), 1, pFile );

	fwrite( &matchingDomainId, sizeof(int), 1, pFile );
	fwrite( &nRotation, sizeof(int), 1, pFile );
	fwrite( &reflectionFlag, sizeof(bool), 1, pFile );
	fwrite( &nHighErrorBin, sizeof(int), 1, pFile );
	if( nHighErrorBin > 0 )
	{
		fwrite( highErrorBinIdList, sizeof(int), nHighErrorBin, pFile );
		fwrite( highErrorList, sizeof(double), nHighErrorBin, pFile );
	}
}

void
FEL_encodedhistogram::read3( FILE* pFile )
{
	fread( &matchingDomainId, sizeof(int), 1, pFile );
	fread( &nRotation, sizeof(int), 1, pFile );
	fread( &reflectionFlag, sizeof(bool), 1, pFile );
}

void
FEL_encodedhistogram::write3( FILE* pFile  )
{
	//fwrite( &spanId, sizeof(int), 1, pFile );

	fwrite( &matchingDomainId, sizeof(int), 1, pFile );
	fwrite( &nRotation, sizeof(int), 1, pFile );
	fwrite( &reflectionFlag, sizeof(bool), 1, pFile );
	//fwrite( &nHighErrorBin, sizeof(int), 1, pFile );
	//if( nHighErrorBin > 0 )
	//{
	//	fwrite( highErrorBinIdList, sizeof(int), nHighErrorBin, pFile );
	//	fwrite( highErrorList, sizeof(double), nHighErrorBin, pFile );
	//}
}


