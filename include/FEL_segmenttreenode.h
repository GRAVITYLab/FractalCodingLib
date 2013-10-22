/*
 * FEL_segmenttreenode.h
 *
 *  Created on: Feb 20, 2012
 *      Author: abon
 */

#ifndef FEL_SEGMENTTREENODE_H_
#define FEL_SEGMENTTREENODE_H_

#include "ITL_header.h"

class FEL_segmenttreenode
{
	int level, nChild;
	int domainID;
	float lowValue, highValue;
	int nMatchedDomains;
	int matchedDomainIDList[50];

	FEL_segmenttreenode *children;
	FEL_segmenttreenode *parent;

public:

	FEL_segmenttreenode();
	FEL_segmenttreenode( float l, float h );

	void initTreeNode( float l, float h, int level );

	void setLevel( int l );
	void setDomainIndex( int i );
	void setRange( float l, float h );
	void getRange( float *bounds );

	void setParent( FEL_segmenttreenode *p );
	void setChildren();

	int getDomainIndex();
	FEL_segmenttreenode* getChild( int index );

	bool isLeaf();

	~FEL_segmenttreenode();
};

#endif
/* FEL_SEGMENTTREENODE_H_ */
