/**
 * Fractal encoding range class
 * Created on: Feb 14, 2012.
 * @author Abon
 * @see Field regular
 */

#ifndef FEL_RANGE_H_
#define FEL_RANGE_H_

#include "ITL_distribution.h"
#include "ITL_field_regular.h"

template <class T>
class FEL_range: public ITL_field_regular<T>
{
	ITL_distribution distribution;

	int rangeID;

	//int timeStamp;
	int nSegment;

	float entropy;

	// For static case
	int bestMatchedDomainId;
	bool ismatchGood;
	double matchError;
	bool isReflectionNeeded;
	int nBinToRotate;

	int nHighErrorBin;
	int* highErrorBinIdList;
	bool* highErrorBinIds;


	// For time varying case
	int* bestMatchedDomainIdList;
	int* bestMatchedSegmentIdList;
	double* matchErrorList;
	int* nBinToRotateList;
	bool* isReflectionNeededList;


public:

	FEL_range()
	{
		//timeStamp = 0;
		nSegment = 1;
		bestMatchedDomainId = -1;
		ismatchGood = true;

		bestMatchedDomainIdList = NULL;
		bestMatchedSegmentIdList = NULL;
		matchErrorList = NULL;
		nBinToRotateList = NULL;
		isReflectionNeededList = NULL;

		bestMatchedDomainIdList = new int[nSegment];
		bestMatchedSegmentIdList = new int[nSegment];
		matchErrorList = new double[nSegment];
		nBinToRotateList = new int[nSegment];
		isReflectionNeededList = new bool[nSegment];
	}

	FEL_range ( int nDim, float *lowSub, float *highSub, int nbin, int ntime = 1 )
	: ITL_field_regular<T>( nDim, lowSub, highSub )
	{
		//timeStamp = 0;
		bestMatchedDomainId = -1;
		nSegment = 1;
		ismatchGood = true;

		//cout << "in range: " << nbin << " " << ntime << endl;
		distribution.initialize( nbin, ntime );

		bestMatchedDomainIdList = NULL;
		bestMatchedSegmentIdList = NULL;
		matchErrorList = NULL;
		nBinToRotateList = NULL;
		isReflectionNeededList = NULL;

		bestMatchedDomainIdList = new int[nSegment];
		bestMatchedSegmentIdList = new int[nSegment];
		matchErrorList = new double[nSegment];
		nBinToRotateList = new int[nSegment];
		isReflectionNeededList = new bool[nSegment];

	}

	FEL_range( const FEL_range<T>& that ) : ITL_field_regular<T>( that )
	{
		this->distribution = that.distribution;
		this->rangeID = that.rangeID;
		this->nSegment = that.nSegment;
		//this->timeStamp = that.timeStamp;
		this->bestMatchedDomainId = that.bestMatchedDomainId;
		this->ismatchGood = that.ismatchGood;
		this->matchError = that.matchError;
		this->isReflectionNeeded = that.isReflectionNeeded;
		this->nBinToRotate = that.nBinToRotate;

		bestMatchedDomainIdList = NULL;
		bestMatchedSegmentIdList = NULL;
		matchErrorList = NULL;
		nBinToRotateList = NULL;
		isReflectionNeededList = NULL;

		//if( this->nSegment > 1 )
		//{
			bestMatchedDomainIdList = new int[nSegment];
			bestMatchedSegmentIdList = new int[nSegment];
			matchErrorList = new double[nSegment];
			nBinToRotateList = new int[nSegment];
			isReflectionNeededList = new bool[nSegment];

			memcpy( this->bestMatchedDomainIdList, that.bestMatchedDomainIdList, sizeof(int)*nSegment );
			memcpy( this->bestMatchedSegmentIdList, that.bestMatchedSegmentIdList, sizeof(int)*nSegment );
			memcpy( this->matchErrorList, that.matchErrorList, sizeof(double)*nSegment );
			memcpy( this->nBinToRotateList, that.nBinToRotateList, sizeof(int)*nSegment );
			memcpy( this->isReflectionNeededList, that.isReflectionNeededList, sizeof(bool)*nSegment );
		//}
	}

	FEL_range( FEL_range<T>& that ) : ITL_field_regular<T>( that )
	{
		this->distribution = that.distribution;
		this->rangeID = that.rangeID;
		this->nSegment = that.nSegment;
		//this->timeStamp = that.timeStamp;
		this->bestMatchedDomainId = that.bestMatchedDomainId;
		this->ismatchGood = that.ismatchGood;
		this->matchError = that.matchError;
		this->isReflectionNeeded = that.isReflectionNeeded;
		this->nBinToRotate = that.nBinToRotate;

		bestMatchedDomainIdList = NULL;
		bestMatchedSegmentIdList = NULL;
		matchErrorList = NULL;
		nBinToRotateList = NULL;
		isReflectionNeededList = NULL;

		//if( this->nSegment > 1 )
		//{
			bestMatchedDomainIdList = new int[nSegment];
			bestMatchedSegmentIdList = new int[nSegment];
			matchErrorList = new double[nSegment];
			nBinToRotateList = new int[nSegment];
			isReflectionNeededList = new bool[nSegment];

			memcpy( this->bestMatchedDomainIdList, that.bestMatchedDomainIdList, sizeof(int)*nSegment );
			memcpy( this->bestMatchedSegmentIdList, that.bestMatchedSegmentIdList, sizeof(int)*nSegment );
			memcpy( this->matchErrorList, that.matchErrorList, sizeof(double)*nSegment );
			memcpy( this->nBinToRotateList, that.nBinToRotateList, sizeof(int)*nSegment );
			memcpy( this->isReflectionNeededList, that.isReflectionNeededList, sizeof(bool)*nSegment );
		//}
	}

	FEL_range& operator= ( const FEL_range<T>& that )
	{
		if (this != &that ) // protect against invalid self-assignment
		{
			this->distribution = that.distribution;
			this->rangeID = that.rangeID;
			this->nSegment = that.nSegment;
			//this->timeStamp = that.timeStamp;
			this->bestMatchedDomainId = that.bestMatchedDomainId;
			this->ismatchGood = that.ismatchGood;
			this->matchError = that.matchError;
			this->isReflectionNeeded = that.isReflectionNeeded;
			this->nBinToRotate = that.nBinToRotate;

			bestMatchedDomainIdList = NULL;
			bestMatchedSegmentIdList = NULL;
			matchErrorList = NULL;
			nBinToRotateList = NULL;
			isReflectionNeededList = NULL;

			//if( this->nSegment > 1 )
			//{
				bestMatchedDomainIdList = new int[nSegment];
				bestMatchedSegmentIdList = new int[nSegment];
				matchErrorList = new double[nSegment];
				nBinToRotateList = new int[nSegment];
				isReflectionNeededList = new bool[nSegment];

				memcpy( this->bestMatchedDomainIdList, that.bestMatchedDomainIdList, sizeof(int)*nSegment );
				memcpy( this->bestMatchedSegmentIdList, that.bestMatchedSegmentIdList, sizeof(int)*nSegment );
				memcpy( this->matchErrorList, that.matchErrorList, sizeof(double)*nSegment );
				memcpy( this->nBinToRotateList, that.nBinToRotateList, sizeof(int)*nSegment );
				memcpy( this->isReflectionNeededList, that.isReflectionNeededList, sizeof(bool)*nSegment );
			//}
		}
		// by convention, always return *this
		return *this;
	}

	FEL_range& operator= ( FEL_range<T>& that )
	{
		if (this != &that ) // protect against invalid self-assignment
		{
			this->distribution = that.distribution;
			this->rangeID = that.rangeID;
			this->nSegment = that.nSegment;
			//this->timeStamp = that.timeStamp;
			this->bestMatchedDomainId = that.bestMatchedDomainId;
			this->ismatchGood = that.ismatchGood;
			this->matchError = that.matchError;
			this->isReflectionNeeded = that.isReflectionNeeded;
			this->nBinToRotate = that.nBinToRotate;

			bestMatchedDomainIdList = NULL;
			bestMatchedSegmentIdList = NULL;
			matchErrorList = NULL;
			nBinToRotateList = NULL;
			isReflectionNeededList = NULL;

			//if( this->nSegment > 1 )
			//{
				bestMatchedDomainIdList = new int[nSegment];
				bestMatchedSegmentIdList = new int[nSegment];
				matchErrorList = new double[nSegment];
				nBinToRotateList = new int[nSegment];
				isReflectionNeededList = new bool[nSegment];

				memcpy( this->bestMatchedDomainIdList, that.bestMatchedDomainIdList, sizeof(int)*nSegment );
				memcpy( this->bestMatchedSegmentIdList, that.bestMatchedSegmentIdList, sizeof(int)*nSegment );
				memcpy( this->matchErrorList, that.matchErrorList, sizeof(double)*nSegment );
				memcpy( this->nBinToRotateList, that.nBinToRotateList, sizeof(int)*nSegment );
				memcpy( this->isReflectionNeededList, that.isReflectionNeededList, sizeof(bool)*nSegment );
			//}
		}
		// by convention, always return *this
		return *this;
	}

	virtual
	~FEL_range()
	{
		if( bestMatchedDomainIdList != NULL ) delete [] bestMatchedDomainIdList;
		if( bestMatchedSegmentIdList != NULL ) delete [] bestMatchedSegmentIdList;
		if( matchErrorList != NULL ) delete [] matchErrorList;
		if( isReflectionNeededList != NULL ) delete [] isReflectionNeededList;
		if( nBinToRotateList != NULL ) delete [] nBinToRotateList;

	}

	void
	setID( int id )
	{
		rangeID = id;
	}

	/*
	void
	setTimeStamp( int t )
	{
		timeStamp = t;
	}
	*/

	void
	setNumSegment( int n )
	{
		nSegment = n;

		fprintf( stderr, "99\n" );
		if( bestMatchedDomainIdList != NULL )	delete [] bestMatchedDomainIdList;
		bestMatchedDomainIdList = new int[nSegment];
		if( bestMatchedSegmentIdList != NULL )	delete [] bestMatchedSegmentIdList;
		bestMatchedSegmentIdList = new int[nSegment];
		if( matchErrorList != NULL )	delete [] matchErrorList;
		matchErrorList = new double[nSegment];
		if( nBinToRotateList != NULL )	delete [] nBinToRotateList;
		nBinToRotateList = new int[nSegment];
		if( isReflectionNeededList != NULL )	delete [] isReflectionNeededList;
		isReflectionNeededList = new bool[nSegment];
		fprintf( stderr, "100\n" );
	}

	void
	setMatchedDomainId( int id )
	{
		bestMatchedDomainId = id;
	}

	void
	setMatQualityFlag( bool f )
	{
		ismatchGood = f;
	}

	void
	setMatchedDomainSegmentId( int* idList, int* sidList )
	{
		assert( idList != NULL );
		//assert( nSegment > 1 );
		memcpy(  bestMatchedDomainIdList, idList, sizeof(int)*nSegment );
		memcpy(  bestMatchedSegmentIdList, sidList, sizeof(int)*nSegment );
	}

	void
	setMatchError( double f )
	{
		matchError = f;
	}

	void
	setMatchError( double* fList )
	{
		assert( fList != NULL );
		//assert( nSegment > 1 );
		memcpy(  matchErrorList, fList, sizeof(double)*nSegment );
	}

	void
	setMatchTransformation( int n, bool f )
	{
		nBinToRotate = n;
		isReflectionNeeded = f;
	}

	void
	setMatchTransformation( int* nList, bool* fList )
	{
		assert( nList != NULL );
		assert( fList != NULL );
		//assert( nSegment > 1 );
		memcpy( nBinToRotateList, nList, sizeof(int)*nSegment );
		memcpy( isReflectionNeededList, fList, sizeof(bool)*nSegment );
	}

	void
	setEntropy( float f )
	{
		entropy = f;
	}

	void
	setDistribution( double *fList, int nbin, int t = 0 )
	{
		assert( fList != NULL );
		distribution.setFrequencies( fList, t );
	}

	void
	setTimevaryingDistribution( double* fList, int nbin, int ntime )
	{
		distribution.setDistribution( fList, nbin, ntime );
	}

	int
	getID()
	{
		return rangeID;
	}

	/*
	int
	getTimeStamp()
	{
		return timeStamp;
	}
	*/

	int
	getNumSegment()
	{
		return nSegment;
	}

	int
	getMatchedDomainId()
	{
		return bestMatchedDomainId;
	}

	bool
	getMatQualityFlag()
	{
		return ismatchGood;
	}

	int
	getMatchedDomainId( int* idList, int* sidList )
	{
		assert( idList != NULL );
		assert( sidList != NULL );
		//assert( nSegment > 1 );
		memcpy( idList, bestMatchedDomainIdList, sizeof(int)*nSegment );
		memcpy( sidList, bestMatchedSegmentIdList, sizeof(int)*nSegment );
	}

	double
	getMatchError()
	{
		return matchError;
	}

	double
	getMatchError( double* eList )
	{
		assert( eList != NULL );
		memcpy( eList, matchErrorList, sizeof(double)*nSegment );
	}

	void
	getMatchTransformation( int* nList, bool* fList )
	{
		assert( nList != NULL );
		assert( fList != NULL );
		//assert( nSegment > 1 );
		memcpy( nList, nBinToRotateList, sizeof(int)*nSegment );
		memcpy( fList, isReflectionNeededList, sizeof(bool)*nSegment );
	}

	float
	getEntropy()
	{
		return entropy;
	}

	void
	getDistribution( double* fList, int* nbin, int t = 0 )
	{
		*nbin = distribution.getNumBin();
		distribution.getFrequencies( fList, t );
	}

	void
	getTimevaryingDistribution( double* fList, int* nbin, int* nt )
	{
		(*nbin) = distribution.getNumBin();
		(*nt) = distribution.getNumTimeStep();
		distribution.getTimeVaryingFrequencies( fList );
	}

	void
	getMatchTransformation( int& n, bool& f )
	{
		n = nBinToRotate;
		f = isReflectionNeeded;
	}

	/**
	 * Range bounds accessor function.
	 */
	//void
	//getBounds( int *limits )
	//{
	//	assert( limits != NULL );
	//	limits[0] = this->grid->lowInt[0];
	//	limits[1] = this->grid->lowInt[1];
	//	limits[2] = this->grid->lowInt[2];
	//	limits[3] = this->grid->highInt[0];
	//	limits[4] = this->grid->highInt[1];
	//	limits[5] = this->grid->highInt[2];
	//}

	static bool
	compare_matcherror( FEL_range<T> ran1, FEL_range<T> ran2 )
	{
		double e1 = ran1.getMatchError();
		double e2 = ran2.getMatchError();

		if( e1 <= e2 ) return true;
		if( e1 > e2 ) return false;

	}// end function

	void
	print()
	{
		int lowInt[4];
		int highInt[4];
		this->getBounds( lowInt, highInt );
		fprintf( stderr, "Extent: <%d %d %d> <%d %d %d> \n", lowInt[0], lowInt[1], lowInt[2],
															 highInt[0], highInt[1], highInt[2] );
		fprintf( stderr, "Matching domain: %d\n", bestMatchedDomainId );
		fprintf( stderr, "Match error: %g\n", matchError );

	}// end function

};

#endif
/* FEL_RANGE_H_ */
