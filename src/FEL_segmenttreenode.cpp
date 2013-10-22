/*
 * FEL_segmenttreenode.cpp
 *
 *  Created on: Feb 20, 2012
 *      Author: abon
 */

#include "FEL_segmenttreenode.h"

FEL_segmenttreenode::FEL_segmenttreenode()
{
	nChild = -1;
	domainID = -1;
	children = NULL;
	parent = NULL;

}// End constructor

FEL_segmenttreenode::FEL_segmenttreenode( float l, float h )
{
	lowValue = l;
	highValue = h;
	nChild = -1;
	domainID = -1;

	children = NULL;
	parent = NULL;

}// End constructor

void
FEL_segmenttreenode::initTreeNode( float l, float h, int lv )
{
	lowValue = l;
	highValue = h;
	nChild = -1;
	domainID = -1;
	level = lv;

	children = NULL;
	parent = NULL;


}// End constructor

void
FEL_segmenttreenode::setLevel( int l )
{
	level = l;

}// End function

void
FEL_segmenttreenode::setDomainIndex( int i )
{
	domainID = i;
}

void
FEL_segmenttreenode::setRange( float l, float h )
{
	lowValue = l;
	highValue = h;

}// End function

void
FEL_segmenttreenode::getRange( float *bounds )
{
	bounds[0] = lowValue;
	bounds[1] = highValue;
}// End function

/**
 *
 */
void
FEL_segmenttreenode::setParent( FEL_segmenttreenode *p )
{
	parent = p;
}// End function

/**
 *
 */
void
FEL_segmenttreenode::setChildren()
{
	nChild = 2;
	children = new FEL_segmenttreenode[2];
}// End function

/**
 *
 */
FEL_segmenttreenode*
FEL_segmenttreenode::getChild( int index )
{
	assert( nChild != 0 );
	return (children + index);
}

/**
 *
 */
int
FEL_segmenttreenode::getDomainIndex()
{
	return domainID;
}

/**
 *
 */
bool
FEL_segmenttreenode::isLeaf()
{
	return (nChild == -1);
}// End function

/**
 *
 */
FEL_segmenttreenode::~FEL_segmenttreenode()
{
	if( children != NULL ) delete [] children;
}// End function


