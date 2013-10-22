/**
 * Fractal encoding domain class
 * Created on: Feb 14, 2012.
 * @author Abon
 * @see Field regular
 */

#ifndef FEL_DOMAIN_CORE_H_
#define FEL_DOMAIN_CORE_H_

#include <cstring>
#include <cassert>
#include <list>

class FEL_domain_core
{
	int nBin;
	double center[3];
	double bounds[6];
	double* freqList;
	float entropy;
	bool isActive;

public:

	FEL_domain_core();
	FEL_domain_core( const FEL_domain_core& that );
	FEL_domain_core( double *limits, int nbin );
	FEL_domain_core( double *limits, double* c, int nbin );
	FEL_domain_core( double *limits, int nbin, double* flist );

	FEL_domain_core& operator= ( const FEL_domain_core& that );

	~FEL_domain_core();

	void getCenter( double* v );
	double getEntropy();
	void getLimits( double* limits );
	double getSize();
	void getDistribution( double* fList );
	bool getActiveStatus();

	void setCenter( double* v );
	void setEntropy( double f );
	void setDistribution( double* fList );
	void setActiveStatus( bool f );

	static bool compare_entropy( FEL_domain_core dom1, FEL_domain_core dom2 );
};

#endif
/* FEL_DOMAIN_CORE_H_ */
