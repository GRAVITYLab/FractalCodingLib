/*
 * FEL_spacetreenode.h
 *
 *  Created on: Feb 23, 2012
 *      Author: abon
 */

#ifndef FEL_SPACETREENODE_H_
#define FEL_SPACETREENODE_H_

#include "ITL_header.h"
#include "FEL_range.h"

template <class T>
class FEL_spacetreenode
{
	int level, nChild;

	float lowLimit[4], highLimit[4];

	FEL_range<T> rangeBlock;
	FEL_spacetreenode<T>* children;
	FEL_spacetreenode<T>* parent;

public:

	FEL_spacetreenode()
	{
		nChild = -1;
		children = NULL;
		parent = NULL;

	}// End constructor

	FEL_spacetreenode( float* l, float* h )
	{
		memcpy( lowLimit, l, sizeof(float)*3 );
		memcpy( highLimit, h, sizeof(float)*3 );

		nChild = -1;

		children = NULL;
		parent = NULL;

	}// End constructor

	void
	initTreeNode( float* l, float* h, int level )
	{
		memcpy( lowLimit, l, sizeof(float)*3 );
		memcpy( highLimit, h, sizeof(float)*3 );

		nChild = -1;
		level = level;

		children = NULL;
		parent = NULL;

	}// End constructor

	void
	setLevel( int l )
	{
		level = l;

	}// End function

	void
	setRangeBlock( int rangeID,
	  			   float* l, float *h,
	  			   double* fList,
	  			   int nbin, int ntime,
	  			   float entropy,
	  			   int rot, bool ref,
	  			   double e, int domID )
	{
		rangeBlock.initialize( 3, l, h );
		rangeBlock.setID( rangeID );
		rangeBlock.setTimevaryingDistribution( fList, nbin, ntime );
		rangeBlock.setEntropy( entropy );
		rangeBlock.setMatchTransformation( rot, ref );
		rangeBlock.setMatchedDomainId( domID );
		rangeBlock.setMatchError( e );
	}

	void
	setRangeMatchQualityFlag( bool flag )
	{
		rangeBlock.setMatQualityFlag( flag );
	}

	void
	setRangeBlockDisrtibution( double* fList, int nbin, int ntime )
	{
		rangeBlock.setTimevaryingDistribution( fList, nbin, ntime );
	}


	void
	setSpatialRange( float* l, float* h )
	{
		memcpy( lowLimit, l, sizeof(float)*3 );
		memcpy( highLimit, h, sizeof(float)*3 );

	}// End function

	void
	getSpatialRange( float* l, float* h )
	{
		memcpy( l, lowLimit, sizeof(float)*3 );
		memcpy( h, highLimit, sizeof(float)*3 );

	}// End function

	void
	getSpatialRange( float *bounds )
	{
		#ifdef DEBUG_MODE
		printf( "Found range: <%g,%g,%g>, <%g,%g,%g>\n", lowLimit[0], lowLimit[1], lowLimit[2],
		 	 	 	 	 	 	 	 	 	 	 	 	highLimit[0], highLimit[1], highLimit[2] );
		#endif

		memcpy( bounds, lowLimit, sizeof(float)*3 );
		memcpy( bounds+3, highLimit, sizeof(float)*3 );

	}// End function

	void setParent( FEL_spacetreenode *p )
	{
		parent = p;

	}// End function

	void setChildren( int nchild )
	{
		nChild = nchild;
		children = new FEL_spacetreenode[nChild];

	}// End function

	int getLevel()
	{
		return level;

	}// End function

	void
	getRangeBlockDisrtibution( double* fList, int nbin, int ntime )
	{
		int nb, nt;
		double* fList2 = new double[nbin*ntime];
		rangeBlock.getTimevaryingDistribution( fList2, &nb, &nt );
		memcpy( fList, fList2, sizeof(double)*(nbin*ntime) );
	}


	FEL_range<T>*
	getRangeBlock()
	{
		return (&rangeBlock);
	}

	FEL_spacetreenode<T>*
	getChild( int index )
	{
		assert( nChild != 0 );
		return (children + index);
	}

	bool
	isLeaf()
	{
		return (nChild == -1);
	}// End function

	void
	print()
	{
		fprintf( stderr, "Level %d: <%g %g %g, %g %g %g>, entropy: %g, match error: %g\n,",
						 level,
						 lowLimit[0], lowLimit[1], lowLimit[2],
						 highLimit[0], highLimit[1], highLimit[2],
						 rangeBlock.getEntropy(),
						 rangeBlock.getMatchError()
						 );
	}

	~FEL_spacetreenode()
	{
		if( children != NULL ) delete [] children;

	}// End function

};



#endif /* FEL_SPACETREENODE_H_ */
