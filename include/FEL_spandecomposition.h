/*
 * FEL_spandecomposition.h
 *
 *  Created on: Aug 10, 2013
 *      Author: abon
 */

#ifndef FEL_SPANDECOMPOSITION_H_
#define FEL_SPANDECOMPOSITION_H_

#include <list>
#include <map>
#include <Spans.h>

class FEL_spandecomposition
{
public:

	static unsigned int unsetRightmostBit( unsigned int x );
	static unsigned int getNextRange_1D( unsigned int x );

	static void decomposeRange_1D( unsigned int x, std::list<Span1D>* partitionList );
	static void decomposeRange_2D( unsigned int x, unsigned int y, std::map<Span2D, int>* xySpanList );
	//static void decomposeRange_3D( unsigned int x, unsigned int y, unsigned int z, std::map<Span3D, int, Span3DComp>* xyzSpanMap );
	static void decomposeRange_3D( unsigned int x, unsigned int y, unsigned int z, std::map<Span3D, int>* xyzSpanMap );
	static void decomposeRange_3D( unsigned int x, unsigned int y, unsigned int z, std::list<Span3D>* xyzSpanList );

};

#endif /* FEL_SPANDECOMPOSITION_H_ */
