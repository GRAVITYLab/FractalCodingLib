#include "FEL_domain_core.h"

FEL_domain_core::FEL_domain_core()
{
	nBin = 0;
	entropy = 0;
	isActive = true;
	freqList = NULL;
}

FEL_domain_core::FEL_domain_core( double *limits, int nbin )
{
	nBin = nbin;
	entropy = 0;
	isActive = true;
	for( int i=0; i<6; i++ )
		bounds[i] = limits[i];

	center[0] = ( bounds[0] + bounds[1] ) / 2.0;
	center[1] = ( bounds[2] + bounds[3] ) / 2.0;
	center[2] = ( bounds[4] + bounds[5] ) / 2.0;

	freqList = new double[nBin];
	memset( freqList, 0, sizeof(double)*nBin );
}

FEL_domain_core::FEL_domain_core( double *limits, double* c, int nbin )
{
	nBin = nbin;
	entropy = 0;
	isActive = true;
	for( int i=0; i<6; i++ )
		bounds[i] = limits[i];

	for( int i=0; i<3; i++ )
		center[i] = c[i];

	freqList = new double[nBin];
	memset( freqList, 0, sizeof(double)*nBin );
}

FEL_domain_core::FEL_domain_core( double *limits, int nbin, double* flist )
{
	nBin = nbin;
	entropy = 0;
	isActive = true;
	for( int i=0; i<6; i++ )
		bounds[i] = limits[i];

	center[0] = ( bounds[0] + bounds[1] ) / 2.0;
	center[1] = ( bounds[2] + bounds[3] ) / 2.0;
	center[2] = ( bounds[4] + bounds[5] ) / 2.0;

	freqList = new double[nBin];
	memcpy( freqList, flist, sizeof(double)*nBin );
}

FEL_domain_core::FEL_domain_core( const FEL_domain_core& that )
{
	this->nBin = that.nBin;
	this->entropy = that.entropy;
	this->isActive = that.isActive;

	this->freqList = NULL;

	for( int i=0; i<6; i++ )
		this->bounds[i] = that.bounds[i];

	for( int i=0; i<3; i++ )
		this->center[i] = that.center[i];


	if( that.freqList != NULL )
	{
		this->freqList = new double[nBin];
		memcpy( this->freqList, that.freqList, sizeof(double)*nBin );
	}
}


FEL_domain_core&
FEL_domain_core::operator= ( const FEL_domain_core& that )
{
	if (this != &that ) // protect against invalid self-assignment
	{
		this->nBin = that.nBin;
		this->entropy = that.entropy;
		this->isActive = that.isActive;

		this->freqList = NULL;
		for( int i=0; i<6; i++ )
			this->bounds[i] = that.bounds[i];

		for( int i=0; i<3; i++ )
			this->center[i] = that.center[i];

		if( that.freqList != NULL )
		{
			this->freqList = new double[nBin];
			memcpy( this->freqList, that.freqList, sizeof(double)*nBin );
		}
	}

	// by convention, always return *this
	return *this;
}

FEL_domain_core::~FEL_domain_core()
{
	if( freqList != NULL )	delete [] freqList;
}

void
FEL_domain_core::getCenter( double* v )
{
	memcpy( v, center, sizeof(double)*3 );
}

double
FEL_domain_core::getEntropy()
{
	return entropy;
}

void
FEL_domain_core::getLimits( double* limits )
{
	for( int i=0; i<6; i++ )
		limits[i] = bounds[i];
}

double
FEL_domain_core::getSize()
{
	return ( ( bounds[1] - bounds[0] + 1 )*
			  ( bounds[3] - bounds[2] + 1 )*
			  ( bounds[5] - bounds[4] + 1 ) );
}


void
FEL_domain_core::getDistribution( double* fList )
{
	assert( fList != NULL );
	memcpy( fList, freqList, sizeof(double)*nBin );
}

bool
FEL_domain_core::getActiveStatus( )
{
	return isActive;
}

void
FEL_domain_core::setCenter( double* v )
{
	memcpy( center, v, sizeof(double)*3 );
}

void
FEL_domain_core::setEntropy( double f )
{
	entropy = f;
}

void
FEL_domain_core::setActiveStatus( bool f )
{
	isActive = f;
}

void
FEL_domain_core::setDistribution( double* fList )
{
	assert( fList != NULL );
	memcpy( freqList, fList, sizeof(double)*nBin );
}

bool
FEL_domain_core::compare_entropy( FEL_domain_core dom1, FEL_domain_core dom2 )
{
	float h1 = dom1.getEntropy();
	float h2 = dom2.getEntropy();

	if( h1 <= h2 ) return true;
	if( h1 > h2 ) return false;

	return false;

}// end function

