/**
 * Fractal encoding domain class
 * Created on: Feb 14, 2012.
 * @author Abon
 * @see Field regular
 */

#ifndef FEL_DOMAIN_H_
#define FEL_DOMAIN_H_

#include <list>

#include "ITL_distribution.h"
#include "ITL_field_regular.h"
#include "FEL_codebookentry.h"

template <class T>
class FEL_domain: public ITL_field_regular<T>
{
	ITL_distribution distribution;

	int nTime;

	float entropy;
	int levelOfSampling;

	//float entropyList[100];
	float* entropyList;

	int nMatchedRanges;
	list<FEL_codebookentry> cbeList;
	//list<int> matchedRangeIDList;
	//list<int> matchedRangeTranslationList;
	//list<bool> matchedRangeReflectionList;
	//list<float> matchedRangeErrorList;

public:

	FEL_domain ()
	{
		nTime = 1;
		nMatchedRanges = 0;
		levelOfSampling = -1;
		entropyList = NULL;
	}

	FEL_domain( const FEL_domain<T>& that ) : ITL_field_regular<T>( that )
	{
		this->nTime = that.nTime;
		this->distribution = that.distribution;
		this->entropy = that.entropy;
		this->levelOfSampling = that.levelOfSampling;
		this->nMatchedRanges = that.nMatchedRanges;
		//this->matchedRangeIDList = that.matchedRangeIDList;
		this->cbeList = that.cbeList;
		if( that.nTime > 1 )
		{
			this->entropyList = new float[nTime];
			memcpy( this->entropyList, that.entropyList, sizeof(float)*nTime );
		}
		else
			entropyList = NULL;
	}

	FEL_domain( FEL_domain<T>& that ) : ITL_field_regular<T>( that )
	{
		this->nTime = that.nTime;
		this->distribution = that.distribution;
		this->entropy = that.entropy;
		this->levelOfSampling = that.levelOfSampling;
		this->nMatchedRanges = that.nMatchedRanges;
		//this->matchedRangeIDList = that.matchedRangeIDList;
		this->cbeList = that.cbeList;
		if( that.nTime > 1 )
		{
			this->entropyList = new float[nTime];
			memcpy( this->entropyList, that.entropyList, sizeof(float)*nTime );
		}
		else
			entropyList = NULL;
	}

	FEL_domain& operator= ( const FEL_domain<T>& that )
	{
		if (this != &that ) // protect against invalid self-assignment
		{
			this->nTime = that.nTime;
			this->distribution = that.distribution;
			this->entropy = that.entropy;
			this->levelOfSampling = that.levelOfSampling;
			this->nMatchedRanges = that.nMatchedRanges;
			//this->matchedRangeIDList = that.matchedRangeIDList;
			this->cbeList = that.cbeList;
			if( that.nTime > 1 )
			{
				this->entropyList = new float[nTime];
				memcpy( this->entropyList, that.entropyList, sizeof(float)*nTime );
			}
			else
				entropyList = NULL;
		}
		// by convention, always return *this
		return *this;
	}

	FEL_domain& operator= ( FEL_domain<T>& that )
	{
		if (this != &that ) // protect against invalid self-assignment
		{
			this->nTime = that.nTime;
			this->distribution = that.distribution;
			this->entropy = that.entropy;
			this->levelOfSampling = that.levelOfSampling;
			this->nMatchedRanges = that.nMatchedRanges;
			//this->matchedRangeIDList = that.matchedRangeIDList;
			this->cbeList = that.cbeList;
			if( that.nTime > 1 )
			{
				this->entropyList = new float[nTime];
				memcpy( this->entropyList, that.entropyList, sizeof(float)*nTime );
			}
			else
				entropyList = NULL;

		}
		// by convention, always return *this
		return *this;
	}

	/**
	 * Copy constructor
	 */
	//FEL_domain ( const FEL_domain<T>& source )
	//{
	//	cout << "here: " << source.nMatchedRanges << endl;
	//	nMatchedRanges = source.nMatchedRanges;
	//	levelOfSampling = source.levelOfSampling;
	//	entropy = source.entropy;
		//matchedRangeIDList = source.matchedRangeIDList;
		//distribution = source.distribution;
	//}

	FEL_domain ( int nDim, float *lowSub, float *highSub, int loD ) : ITL_field_regular<T>( nDim, lowSub, highSub )
	{
		nMatchedRanges = 0;
		levelOfSampling = loD;

		nTime = 1;
		entropyList = NULL;
	}

	FEL_domain ( int nDim, float *lowSub, float *highSub, int loD, int nbin, int ntime = 1 )
	: ITL_field_regular<T>( nDim, lowSub, highSub )
	{
		nMatchedRanges = 0;
		levelOfSampling = loD;
		nTime = ntime;

		distribution.initialize( nbin, ntime );

		if( nTime > 1 )
			entropyList = new float[nTime];
		else
			entropyList = NULL;
	}

	virtual
	~FEL_domain()
	{
		if( entropyList != NULL )	delete [] entropyList;
	}

	/**
	 * Stores matching range information
	 */
	/*
	void
	addMatchingRange( int rangeID )
	{
		//matchedRangeIDList.push_back( rangeID );

		nMatchedRanges ++;
	}

	*/

	void
	addMatchingRange( int rangeID, float matchError, int nBinToRot, bool isRef )
	{
		//matchedRangeIDList.push_back( rangeID );
		//matchedRangeTranslationList.push_back( nBinToRotate );
		//matchedRangeReflectionList.push_back( isRef );
		//matchedRangeErrorList.push_back( matchError );

		FEL_codebookentry newEntry( rangeID, matchError, nBinToRot, isRef );
		cbeList.push_back( newEntry );

		nMatchedRanges ++;
	}

	void
	addMatchingRangeSegment( int rangeID, int rangeSegmentID,
							 int bestMatchedSegmentID,
							 float matchError,
							 int nBinToRot,
							 bool isRef )
	{
		FEL_codebookentry newEntry( rangeID, rangeSegmentID, bestMatchedSegmentID, matchError, nBinToRot, isRef );
		cbeList.push_back( newEntry );

		nMatchedRanges ++;
	}

	void
	sortCodeBookEntries()
	{
		cbeList.sort( FEL_codebookentry::compare_matcherror );
	}


	/**
	 * Returns underlying data level
	 */
	int
	getDataLoD()
	{
		return levelOfSampling;
	}

	/**
	 * Returns number of matching range information
	 */
	int
	getNumMatchingRange()
	{
		return nMatchedRanges;
	}

	/**
	 * Returns matching range information
	 */
	void
	getMatchingRange( int *rangeList )
	{
		int i = 0;
		//for( list<int>::iterator rIter = matchedRangeIDList.begin();
		//	 rIter != matchedRangeIDList.end(); ++rIter )
		//{
		//	rangeList[i] = (*rIter);
		//	i++;
		//}
		for( list<FEL_codebookentry>::iterator cbIter = cbeList.begin();
			 cbIter != cbeList.end(); ++cbIter )
		{
			rangeList[i] = cbIter->getRangeID();
			i++;
		}


	}

	int
	getMatchingRangeInformation( int iD, double* error, int* nRot, bool* isRef )
	{
		//list<int>::iterator rIter = matchedRangeIDList.begin();
		//list<int>::iterator t1Iter = matchedRangeTranslationList.begin();
		//list<bool>::iterator t2Iter = matchedRangeReflectionList.begin();
		list<FEL_codebookentry>::iterator cbIter = cbeList.begin();

		for( int i=0; i<iD; i++ )
			cbIter ++;

		(*nRot) = cbIter->getRotation();
		(*isRef) = cbIter->getReflection();
		(*error) = cbIter->getMatchError();

		return cbIter->getRangeID();
	}

	/**
	 * Domain bounds accessor function.
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

	//int
	//getSize()
	//{
	//	return this->grid->nVertices;
	//}

	float
	getEntropy()
	{
		return entropy;
	}

	void
	getTimevaryingEntropy( float* hList )
	{
		assert( hList != NULL );
		//assert( entropyList != NULL );
		memcpy( hList, entropyList, sizeof(float)*nTime );

	}

	void
	getDistribution( double* fList, int* nbin, int t = 0 )
	{
		*nbin = distribution.getNumBin();
		distribution.getFrequencies( fList, t );
	}

	void
	getTimeVaryingDistribution( double* fList, int* nbin, int* nt )
	{
		(*nbin) = distribution.getNumBin();
		(*nt) = distribution.getNumTimeStep();
		distribution.getTimeVaryingFrequencies( fList );
	}

	void
	getTimeVaryingDistributionOfSegment( int segStart, int segLength, double* fList, int* nbin, int* nt )
	{
		(*nbin) = distribution.getNumBin();
		(*nt) = distribution.getNumTimeStep();

		double temp[(*nbin)*(*nt)];
		distribution.getTimeVaryingFrequencies( temp );

		memcpy( fList, temp+(*nbin)*(segStart), sizeof(double)*((*nbin)*segLength) );

	}

	void
	setEntropy( float f )
	{
		entropy = f;
	}

	void
	setTimevaryingEntropy( float* hList )
	{
		assert( hList != NULL );
		//assert( entropyList != NULL );
		cout << nTime << endl;
		memcpy( entropyList, hList, sizeof(float)*nTime );
	}


	void
	setDistribution( double* fList, int nbin, int t = 0 )
	{
		assert( fList != NULL );
		distribution.setDistribution( fList, nbin );
	}

	void
	setTimevaryingDistribution( double* fList, int nbin, int time )
	{
		distribution.setDistribution( fList, nbin, time );
	}

	static bool
	compare_entropy( FEL_domain<T> dom1, FEL_domain<T> dom2 )
	{
		float h1 = dom1.getEntropy();
		float h2 = dom2.getEntropy();

		if( h1 <= h2 ) return true;
		if( h1 > h2 ) return false;

	}// end function

};
#endif
/* FEL_DOMAIN_H_ */
