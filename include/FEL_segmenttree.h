/*
 * FEL_segmenttree.h
 *
 *  Created on: Feb 20, 2012
 *      Author: abon
 */

#ifndef FEL_SEGMENTTREE_H_
#define FEL_SEGMENTTREE_H_

#include <list>
#include "FEL_domain.h"
#include "FEL_segmenttreenode.h"

#define EPSILON1 0.01

template <class T>
class FEL_segmenttree
{
	int nChild;

	//  The list of all domains, needed for tree construction.
	list< FEL_domain<T> > domainList;
	int nDomainTotal;

	FEL_segmenttreenode *root;
	FEL_segmenttreenode *curNode;

public:

	/**
	 * Constructor
	 */
	FEL_segmenttree()
	{
		root = NULL;
		curNode = NULL;
	}

	/**
	 * Constructor
	 */
	FEL_segmenttree( list< FEL_domain<T> > dList )
	{
		domainList = dList;
		nDomainTotal = domainList.size();

		root = NULL;
		curNode = NULL;
	}

	/**
	 * Initialization of tree.
	 * Create a single node tree that contains the entire range
	 */
	void
	initTree( list< FEL_domain<T> > dList )
	{
		// Get domain list to place in the tree
		domainList = dList;

		// Create root node using the first and last nodes
		// of the domain List
		float low = getListItem( 0 );
		float high = getListItem( nDomainTotal-1 );
		root = new FEL_segmenttreenode( low, high );
		curNode = root;

		// Set tree properties
		nChild = 2;

	}// End constructor

	/**
	 * Recursive function that starts from the root and
	 * subdivides the value range into segments.
	 */
	void createTree( FEL_segmenttreenode *cur,
					 int lowIndex, int highIndex,
					 int level )
	{
		// Add range to the node
		// Terminate iteration if the two ends of
		// the range has converged
		if( lowIndex >= highIndex )
		{
			if( lowIndex > highIndex )	lowIndex = highIndex;

			float v = getListItem( lowIndex );
			cur->initTreeNode( v, v, level );
			cur->setDomainIndex( lowIndex );

			#ifdef DEBUG_MODE
			float lim[2];
			cur->getRange( lim );
			fprintf( stderr, "Segment tree node: %d: <%d: %g, %d: %g>\n",
							 level, lowIndex, lim[0], highIndex, lim[1] );
			#endif

			return;
		}
		else
		{
			float lv = getListItem( lowIndex );
			float hv = getListItem( highIndex );
			cur->initTreeNode( lv, hv, level );
		}

		#ifdef DEBUG_MODE
		float lim[2];
		cur->getRange( lim );
		fprintf( stderr, "Segment tree Leaf node: %d: <%d: %g, %d: %g>\n",
						 level, lowIndex, lim[0], highIndex, lim[1] );
		#endif

		// Further subdivide the range
		int mid = (lowIndex+highIndex) / 2;

		// Add children to the current node
		cur->setChildren();

		// Recursive call to children
		createTree( cur->getChild(0), lowIndex, mid, level+1 );
		createTree( cur->getChild(1), mid+1, highIndex, level+1 );

	}// End function

	list<int>
	searchTree( FEL_segmenttreenode *cur,
				float qLow, float qHigh )
	{
		list<int> localSearchResultList;

		assert( qLow <= qHigh );

		// Get values associated to this node
		float lim[2];
		cur->getRange( lim );
		//printf( "At Node: <%g,%g>\n", lim[0], lim[1] );

		// Return if arrived at leaf node
		if( cur->isLeaf() )
		{
			//printf( "At leaf: <%g,%g>\n", lim[0], lim[1] );

			// Print result if it is winthin
			// query range
			if( qLow <= lim[0] && lim[1] <= qHigh )
			{
				int domainID = cur->getDomainIndex();
				printf( "Returned entropy and domain ID: <%g,%g>, %d\n", lim[0], lim[1], domainID );
				localSearchResultList.push_back( domainID );
			}

			return localSearchResultList;
		}

		// If the node's range and the query range has
		// no overlap, then stop recursing
		if( qHigh < lim[0] || qLow > lim[1] )
		{
			//printf( "not recursing down: <%g,%g>\n", lim[0], lim[1] );
			return localSearchResultList;
		}
		// If the node's range is contained by the query range,
		// If the node's range contains the query range,
		// or they have a partial overlap,
		// need to recurse down to finer nodes
		else if( ( qLow <= lim[0] && lim[1] <= qHigh ) ||
				 ( lim[0] <= qLow && qHigh <= lim[1] ) ||
				 ( qLow <= lim[0] && qHigh >= lim[0] && qHigh <= lim[1] ) ||
				 ( lim[0] <= qLow && qLow <= lim[1] && lim[1] <= qHigh ) )
		{
			//printf( "recursing down: <%g,%g>\n", lim[0], lim[1] );
			list<int> c1List = searchTree( cur->getChild(0), qLow, qHigh );
			list<int> c2List = searchTree( cur->getChild(1), qLow, qHigh );

			list<int>::iterator it = localSearchResultList.end();
			if( !c1List.empty() )
				localSearchResultList.insert( it, c1List.begin(), c1List.end() );
			it = localSearchResultList.end();
			if( !c2List.empty() )
				localSearchResultList.insert( it, c2List.begin(), c2List.end() );
		}

		return localSearchResultList;

	}// End function

	FEL_segmenttreenode*
	getRoot()
	{
		return root;
	}

	float
	getListItem( int index )
	{
		int i=0;
		for( typename list< FEL_domain<T> >::iterator domainIter = domainList.begin();
			 domainIter != domainList.end();
			 domainIter++ )
		{
			if( i == index )
				return domainIter->getEntropy();
			else
				i++;

		}// end for

		return -1.0f;

	}// End function

};


#endif /* FEL_SEGMENTTREE_H_ */
