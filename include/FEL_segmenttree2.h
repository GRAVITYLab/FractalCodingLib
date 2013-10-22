/*
 * FEL_segmenttree2.h
 *
 *  Created on: Sep 15, 2012
 *      Author: abon
 */

#ifndef FEL_SEGMENTTREE2__H_
#define FEL_SEGMENTTREE2_H_

#include <list>
#include "FEL_domain_core.h"
#include "FEL_segmenttreenode.h"

#define EPSILON1 0.01

class FEL_segmenttree2
{
	int nChild;

	//  The list of all domains, needed for tree construction.
	list< FEL_domain_core> domainList;
	int nDomainTotal;

	FEL_segmenttreenode *root;
	FEL_segmenttreenode *curNode;

public:

	FEL_segmenttree2();
	FEL_segmenttree2( list< FEL_domain_core> dList );

	void initTree( list< FEL_domain_core> dList );
	void createTree( FEL_segmenttreenode *cur, int lowIndex, int highIndex, int level );
	list<int> searchTree( FEL_segmenttreenode *cur, float qLow, float qHigh );
	FEL_segmenttreenode* getRoot();
	float getListItem( int index );
};

#endif

/* FEL_SEGMENTTREE2_H_ */
