/*
 * FEL_encodedSimpleHistogram.h
 *
 *  Created on: Mar 10, 2013
 *      Author: abon
 */

#ifndef FEL_ENCODEDSIMPLEHISTOGRAM_H_
#define FEL_ENCODEDSIMPLEHISTOGRAM_H_

#include <cstdio>
#include <cassert>
#include <string>
#include <cstring>
#include <list>

using namespace std;

class FEL_encodedsimplehistogram
{
	int low[3], high[3];
	int spanId;
	int nNonZeroBin;
	int* binIdList;
	double* freqList;

public:

	FEL_encodedsimplehistogram();
	FEL_encodedsimplehistogram( 	double* fList, int nBin );
	FEL_encodedsimplehistogram( 	double* fList, int* l, int* h, int nBin );

	FEL_encodedsimplehistogram& operator= ( const FEL_encodedsimplehistogram& that );
	FEL_encodedsimplehistogram( const FEL_encodedsimplehistogram& that );

	//Destructor
	~FEL_encodedsimplehistogram();

	int getNumNonZeroBins();
	void getSpan( int* l, int* h );
	int getSpan2( int* l, int* h );

	void setSpan( int* l, int* h );
	void setSpan( int* l, int* h, int id );

	void encode( double* fList, int nBin );
	void decode( double* fList, int nBin );

	void read( FILE* pFile );
	void readSpan( FILE* pFile );
	void readCount( FILE* pFile );
	void readBinIds( FILE* pFile );
	void readBinFrequencies( FILE* pFile );

	void print();
	void write( FILE* pFile );
	void writeSpan( FILE* pFile );
	void writeCount( FILE* pFile );
	void writeBinIds( FILE* pFile );
	void writeBinFrequencies( FILE* pFile );
};

#endif /* FEL_ENCODEDSIMPLEHISTOGRAM_H_ */
