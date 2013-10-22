/*
 * FEL_spandecomposition.cpp
 *
 *  Created on: Aug 10, 2013
 *      Author: abon
 */

#include "FEL_spandecomposition.h"

unsigned int
FEL_spandecomposition::unsetRightmostBit( unsigned int x )
{
	x = x - (x&-x);
	return x;
}

unsigned int
FEL_spandecomposition::getNextRange_1D( unsigned int x )
{
	unsigned int xLow = unsetRightmostBit( x );

	#if DEBUG_MODE
	fprintf( stderr, "%u becomes %u after bit unsetting\n", x, xLow );
	#endif

	return ( x - xLow );
}

void
FEL_spandecomposition::decomposeRange_1D( unsigned int x, std::list<Span1D>* spanList )
{
	int input = x;
	unsigned int xLow = 0, residue = 0;

	do
	{
		xLow = getNextRange_1D( x );

		residue = x - xLow;

		Span1D nextRange( residue+1, x );
		spanList->push_back( nextRange );

		x = residue;

	}while( residue > 1 );

	#ifdef DEBUG_MODE
	fprintf( stderr, "Decomposition for %d:\n", input );
	for( std::list<Span1D>::iterator iter = spanList->begin();
		  iter != spanList->end();
		  iter++ )
		fprintf( stderr, "[%u %u]\n", iter->low, iter->high );
	#endif
}

void
FEL_spandecomposition::decomposeRange_2D( unsigned int x, unsigned int y,
												std::map<Span2D, int>* xySpanMap )
{
	std::list<Span1D> xPartitionList;
	std::list<Span1D> yPartitionList;

	decomposeRange_1D( x, &xPartitionList );
	decomposeRange_1D( y, &yPartitionList );

	int nXPart = xPartitionList.size();
	int nYPart = yPartitionList.size();

	std::list<Span1D>::iterator iterx = xPartitionList.begin();
	for( int x=0; x<nXPart; x++ )
	{
		std::list<Span1D>::iterator itery = yPartitionList.begin();
		for( int y=0; y<nYPart; y++ )
		{
			Span2D nextRange( iterx->low, iterx->high,
							   itery->low, itery->high );

			xySpanMap->insert( std::pair<Span2D,int>( nextRange, 0 ) );

			itery++;
		}
		iterx++;
	}

	#ifdef DEBUG_MODE
	for( std::map<Span2D, int>::iterator iter = xySpanMap->begin();
		  iter != xyMap->end();
		  iter++ )
		fprintf( stderr, "[%u %u]-[%u %u]\n", iter->first.low[0], iter->first.high[0],
											  iter->first.low[1], iter->first.high[1] );
	#endif

	xPartitionList.clear();
	yPartitionList.clear();
}

void
FEL_spandecomposition::decomposeRange_3D( unsigned int x, unsigned int y, unsigned int z,
												//std::map<Span3D, int, Span3DComp>* xyzSpanMap )
		std::map<Span3D, int>* xyzSpanMap )
{
	std::list<Span1D> xPartitionList;
	std::list<Span1D> yPartitionList;
	std::list<Span1D> zPartitionList;

	decomposeRange_1D( x, &xPartitionList );
	decomposeRange_1D( y, &yPartitionList );
	decomposeRange_1D( z, &zPartitionList );

	int nXPart = xPartitionList.size();
	int nYPart = yPartitionList.size();
	int nZPart = zPartitionList.size();
	#if DEBUG_MODE
	fprintf( stderr, "[%d, %d, %d]\n", x, y, z );
	fprintf( stderr, "1D lists contains [%d, %d, %d] spans\n", nXPart, nYPart, nZPart );
	#endif

	std::list<Span1D>::iterator iterx = xPartitionList.begin();
	for( int x=0; x<nXPart; x++ )
	{
		std::list<Span1D>::iterator itery = yPartitionList.begin();
		for( int y=0; y<nYPart; y++ )
		{
			std::list<Span1D>::iterator iterz = zPartitionList.begin();
			for( int z=0; z<nZPart; z++ )
			{
				Span3D nextRange( iterx->low, iterx->high,
								   itery->low, itery->high,
								   iterz->low, iterz->high );

				std::pair<std::map<Span3D,int>::iterator,bool> ret = xyzSpanMap->insert( std::pair<Span3D,int>( nextRange, 0) );
				#if DEBUG_MODE
				if( ret.second == false )
				{
					nextRange.print();
					fprintf( stderr, "Already exists\n" );
				}
				else
				{
					nextRange.print();
					fprintf( stderr, "inserted.\n" );
				}
				#endif

				iterz++;
			}
			itery++;
		}
		iterx++;
	}

	#ifdef DEBUG_MODE
	for( std::map<Span3D, int>::iterator iter = xyzSpanMap->begin();
		  iter != xyzSpanMap->end();
		  iter++ )
		fprintf( stderr, "[%u %u]-[%u %u]-[%u %u]\n", iter->first.low[0], iter->first.high[0],
											   	   	   iter->first.low[1], iter->first.high[1],
											   	   	   iter->first.low[2], iter->first.high[2] );
	#endif

	xPartitionList.clear();
	yPartitionList.clear();
	zPartitionList.clear();
}

void
FEL_spandecomposition::decomposeRange_3D( unsigned int x, unsigned int y, unsigned int z,
												std::list<Span3D>* xyzSpanList )
{
	std::list<Span1D> xPartitionList;
	std::list<Span1D> yPartitionList;
	std::list<Span1D> zPartitionList;

	decomposeRange_1D( x, &xPartitionList );
	decomposeRange_1D( y, &yPartitionList );
	decomposeRange_1D( z, &zPartitionList );

	int nXPart = xPartitionList.size();
	int nYPart = yPartitionList.size();
	int nZPart = zPartitionList.size();
	#ifdef DEBUG_MODE
	fprintf( stderr, "1D lists contains [%d, %d, %d] spans\n", nXPart, nYPart, nZPart );
	#endif

	std::list<Span1D>::iterator iterx = xPartitionList.begin();
	for( int x=0; x<nXPart; x++ )
	{
		std::list<Span1D>::iterator itery = yPartitionList.begin();
		for( int y=0; y<nYPart; y++ )
		{
			std::list<Span1D>::iterator iterz = zPartitionList.begin();
			for( int z=0; z<nZPart; z++ )
			{
				Span3D nextRange( iterx->low, iterx->high,
								   itery->low, itery->high,
								   iterz->low, iterz->high );

				xyzSpanList->push_back( nextRange );

				iterz++;
			}
			itery++;
		}
		iterx++;
	}

	xPartitionList.clear();
	yPartitionList.clear();
	zPartitionList.clear();
}


/*
void
decomposeRange_3D_2( int nX, int nY, int nZ )
{
	for( int z = 0; z<=nX-1; z++ )
	{
		for( int y = 0; y<=nY-1; y++ )
		{
			for( int x = 0; x<=nX-1; x++ )
			{
				doneField[x][y][z] = false;
			}
		}
	}

	for( int z = nZ; z>=1; z-- )
	{
		for( int y = nY; y>=1; y-- )
		{
			for( int x = nX; x>=1; x-- )
			{
				if( !doneField[x-1][y-1][z-1] )
					decomposeRange_3D_rec( x, y, z );
			}
		}
	}

}// end function

// 1 to N range
void
decomposeRange_3D_rec( int x, int y, int z )
{
	if( x == 0 && y == 0 && z == 0 )
		return;

	if( doneField[x][y][z] == true )
		return;

	int xPart = 0, yPart = 0, zPart = 0;

	// Find the latest span for this location
	xPart = getNextRange_1D( x );
	yPart = getNextRange_1D( y );
	zPart = getNextRange_1D( z );

	Span3D nextSpan( x - xPart + 1, x,
					 y - yPart + 1, y,
					 z - zPart + 1, z );
	xyzSpanList.push_back( nextSpan );

	//MapValue nextMap( nextSpan.getSize(), 0 );
	//std::pair<Span3D,MapValue> nextEntry( nextSpan, nextMap );
	//part3dMap.insert( nextEntry );

	doneField[x][y][z] = true;

	decomposeRange_3D_2( x-xPart, y, z );
	decomposeRange_3D_2( x, y-yPart, z );
	decomposeRange_3D_2( x, y, z-zPart );

	decomposeRange_3D_2( x-xPart, y-yPart, z );
	decomposeRange_3D_2( x, y-yPart, z-zPart );
	decomposeRange_3D_2( x-xPart, y, z-zPart );

	decomposeRange_3D_2( x-xPart, y-yPart, z-zPart );

}// end function
*/



