/*
 * FEL_encodedHistogram.h
 *
 *  Created on: Sep 9, 2012
 *      Author: abon
 */

#ifndef FEL_ENCODEDHISTOGRAM_H_
#define FEL_ENCODEDHISTOGRAM_H_

#include <cstdio>
#include <cassert>
#include <string>
#include <cstring>
#include <list>

using namespace std;

class FEL_encodedhistogram
{
	int low[3];
	int high[3];
	int spanId;

	int matchingDomainId;
	int nRotation;
	bool reflectionFlag;
	long nElement;
	double encodingError;
	int nHighErrorBin;
	int* highErrorBinIdList;
	double* highErrorList;

public:

	FEL_encodedhistogram();
	FEL_encodedhistogram(	 int domainId, int nElem, int nHighErrorBin,
							 int* highErrorBinId, double* highErrorList );
	FEL_encodedhistogram( int* low, int* high,
							 int domainId, int nElem, int nHighErrorBin,
							 int* highErrorBinId, double* highErrorList );

	FEL_encodedhistogram& operator= ( const FEL_encodedhistogram& that );
	FEL_encodedhistogram( const FEL_encodedhistogram& that );

	//Destructor
	~FEL_encodedhistogram();

	void getSpan( int* l, int* h );
	int getSpan2( int* l, int* h );
	int getSpanId();
	int getMatchingDomainId();
	long getNumDataElement();
	double getEncodingError();
	int getNumHighErrorBin();
	int getRotation();
	bool getReflection();
	void getHighErrorBinInfo( int* binIds, double* highErrors );
	void getHighErrorBinIds( int* binIds );
	void getHighErrorBinFreqs( double* highErrors );

	void setSpan( int* l, int* h );
	void setSpan( int* l, int* h, int id );
	void setSpanId( int id );
	void setMatchingDomainId( int n );
	void setNumDataElement( long n );
	void setEncodingError( double d );
	void setNumHighErrorBin( int n );
	void setRotation( int nrot );
	void setReflection( bool isref );
	void setHighErrorBinInfo( int* binIds, double* highErrors );
	void setHighErrorBinIds( int* binIds );
	void setHighErrorBinFreqs( double* highErrors );

	void read( FILE* pFile );
	void write( FILE* pFile );

	void read2( FILE* pFile );
	void write2( FILE* pFile );

	void read3( FILE* pFile );
	void write3( FILE* pFile );
};

#endif /* FEL_encodedhistogram_H_ */
