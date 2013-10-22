/*
 * FEL_encodedSimpleHistogram.cpp
 *
 *  Created on: Mar 10, 2013
 *      Author: abon
 */

#include "FEL_encodedsimplehistogram.h"

FEL_encodedsimplehistogram::FEL_encodedsimplehistogram()
{
	// Zero-set range
	memset( low, 0, sizeof(int)*3 );
	memset( high, 0, sizeof(int)*3 );

	nNonZeroBin = 0;
	spanId = 0;

	binIdList = 0;
	freqList = 0;
}

FEL_encodedsimplehistogram::FEL_encodedsimplehistogram( double* fList, int nBin )
{
	// Zero-set range
	memset( low, 0, sizeof(int)*3 );
	memset( high, 0, sizeof(int)*3 );

	// Count the non-zero bins
	nNonZeroBin = 0;
	for( int i=0; i<nBin; i++ )
	{
		if( fList[i] != 0 )
			nNonZeroBin = 0;
	}
	spanId = 0;

	// Set the non-zero bins
	binIdList = new int[nNonZeroBin];
	freqList = new double[nNonZeroBin];
	int id = 0;
	for( int i=0; i<nBin; i++ )
	{
		if( fList[i] != 0 )
		{
			binIdList[id] = i;
			freqList[id] = fList[i];
			id++;
		}
	}

}// End constructor

FEL_encodedsimplehistogram::FEL_encodedsimplehistogram( double* fList,
															   int* l,
															   int* h,
															   int nBin )
{
	// Set range
	memcpy( low, l, sizeof(int)*3 );
	memcpy( high, h, sizeof(int)*3 );

	// Count the non-zero bins
	nNonZeroBin = 0;
	for( int i=0; i<nBin; i++ )
	{
		if( fList[i] != 0 )
			nNonZeroBin = 0;
	}
	spanId = 0;

	// Set the non-zero bins
	binIdList = new int[nNonZeroBin];
	freqList = new double[nNonZeroBin];
	int id = 0;
	for( int i=0; i<nBin; i++ )
	{
		if( fList[i] != 0 )
		{
			binIdList[id] = i;
			freqList[id] = fList[i];
			id++;
		}
	}

}// End constructor


FEL_encodedsimplehistogram::FEL_encodedsimplehistogram( const FEL_encodedsimplehistogram& that )
{
	// Set range
	memcpy( this->low, that.low, sizeof(int)*3 );
	memcpy( this->high, that.high, sizeof(int)*3 );

	this->nNonZeroBin = that.nNonZeroBin;
	this->spanId = that.spanId;

	this->binIdList = 0;
	this->freqList = 0;

	if( this->nNonZeroBin > 0 )
	{
		this->binIdList = new int[nNonZeroBin];
		this->freqList = new double[nNonZeroBin];

		memcpy( this->binIdList, that.binIdList, sizeof(int)*nNonZeroBin );
		memcpy( this->freqList, that.freqList, sizeof(double)*nNonZeroBin );
	}

}

FEL_encodedsimplehistogram&
FEL_encodedsimplehistogram::operator= ( const FEL_encodedsimplehistogram& that )
{
	if (this != &that ) // protect against invalid self-assignment
	{
		// Set range
		memcpy( this->low, that.low, sizeof(int)*3 );
		memcpy( this->high, that.high, sizeof(int)*3 );

		this->nNonZeroBin = that.nNonZeroBin;
		this->spanId = that.spanId;

		this->binIdList = 0;
		this->freqList = 0;

		if( this->nNonZeroBin > 0 )
		{
			this->binIdList = new int[nNonZeroBin];
			this->freqList = new double[nNonZeroBin];

			memcpy( this->binIdList, that.binIdList, sizeof(int)*nNonZeroBin );
			memcpy( this->freqList, that.freqList, sizeof(double)*nNonZeroBin );
		}
	}

	// by convention, always return *this
	return *this;
}

FEL_encodedsimplehistogram::~FEL_encodedsimplehistogram()
{
	if( binIdList != 0 )	delete [] binIdList;
	if( freqList != 0 )	delete [] freqList;
}

int
FEL_encodedsimplehistogram::getNumNonZeroBins()
{
	return nNonZeroBin;
}

void
FEL_encodedsimplehistogram::getSpan( int* l, int* h )
{
	memcpy( l, low, sizeof(int)*3 );
	memcpy( h, high, sizeof(int)*3 );
}

int
FEL_encodedsimplehistogram::getSpan2( int* l, int* h )
{
	memcpy( l, low, sizeof(int)*3 );
	memcpy( h, high, sizeof(int)*3 );
	return spanId;
}

void
FEL_encodedsimplehistogram::setSpan( int* l, int* h )
{
	memcpy( low, l, sizeof(int)*3 );
	memcpy( high, h, sizeof(int)*3 );
}

void
FEL_encodedsimplehistogram::setSpan( int* l, int* h, int id )
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
FEL_encodedsimplehistogram::encode( double* fList, int nBin )
{
	// Count the non-zero bins
	nNonZeroBin = 0;
	for( int i=0; i<nBin; i++ )
	{
		if( fList[i] != 0 )
			nNonZeroBin ++;
	}

	// Set the non-zero bins
	binIdList = new int[nNonZeroBin];
	freqList = new double[nNonZeroBin];
	int id = 0;
	for( int i=0; i<nBin; i++ )
	{
		if( fList[i] != 0 )
		{
			binIdList[id] = i;
			freqList[id] = fList[i];
			id++;
		}
	}

}

void
FEL_encodedsimplehistogram::decode( double* fList, int nBin )
{
	memset( fList, 0.0, sizeof(double)*nBin );

	// Set the non-zero bins
	for( int i=0; i<nNonZeroBin; i++ )
	{
		fList[binIdList[i]] = freqList[i];
	}
}

void
FEL_encodedsimplehistogram::print()
{
	for( int i=0; i<nNonZeroBin; i++ )
	{
		fprintf( stderr, "[%d: %g], ", binIdList[i], freqList[i] );
	}
	fprintf( stderr, "\n" );
}

void
FEL_encodedsimplehistogram::read( FILE* pFile )
{
	fread( low, sizeof(int), 3, pFile );
	fread( high, sizeof(int), 3, pFile );

	fread( &nNonZeroBin, sizeof(int), 1, pFile );

	if( nNonZeroBin > 0 )
	{
		fread( binIdList, sizeof(int), nNonZeroBin, pFile );
		fread( freqList, sizeof(double), nNonZeroBin, pFile );
	}
}

void
FEL_encodedsimplehistogram::readSpan( FILE* pFile )
{
	fread( low, sizeof(int), 3, pFile );
	fread( high, sizeof(int), 3, pFile );
}

void
FEL_encodedsimplehistogram::readCount( FILE* pFile )
{
	fread( &nNonZeroBin, sizeof(int), 1, pFile );
}

void
FEL_encodedsimplehistogram::readBinIds( FILE* pFile )
{
	if( nNonZeroBin > 0 )
	{
		if( binIdList == 0 )
			binIdList = new int[nNonZeroBin];

		fread( binIdList, sizeof(int), nNonZeroBin, pFile );
	}
}

void
FEL_encodedsimplehistogram::readBinFrequencies( FILE* pFile )
{
	if( nNonZeroBin > 0 )
	{
		if( freqList == 0 )
			freqList = new double[nNonZeroBin];

		fread( freqList, sizeof(double), nNonZeroBin, pFile );
	}
}

void
FEL_encodedsimplehistogram::write( FILE* pFile )
{
	fwrite( low, sizeof(int), 3, pFile );
	fwrite( high, sizeof(int), 3, pFile );

	fwrite( &nNonZeroBin, sizeof(int), 1, pFile );

	if( nNonZeroBin > 0 )
	{
		fwrite( binIdList, sizeof(int), nNonZeroBin, pFile );
		fwrite( freqList, sizeof(double), nNonZeroBin, pFile );
	}
}

void
FEL_encodedsimplehistogram::writeSpan( FILE* pFile )
{
	fwrite( low, sizeof(int), 3, pFile );
	fwrite( high, sizeof(int), 3, pFile );

}

void
FEL_encodedsimplehistogram::writeCount( FILE* pFile )
{
	fwrite( &nNonZeroBin, sizeof(int), 1, pFile );
}

void
FEL_encodedsimplehistogram::writeBinIds( FILE* pFile )
{
	if( nNonZeroBin > 0 )
		fwrite( binIdList, sizeof(int), nNonZeroBin, pFile );
}

void
FEL_encodedsimplehistogram::writeBinFrequencies( FILE* pFile )
{
	if( nNonZeroBin > 0 )
		fwrite( freqList, sizeof(double), nNonZeroBin, pFile );
}





