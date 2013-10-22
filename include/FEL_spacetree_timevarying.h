/*
 * FEL_spacetree_timevarying.h
 *
 *  Created on: Mar 18, 2012
 *      Author: abon
 */

#ifndef FEL_SPACETREE_TIMEVARYING_H_
#define FEL_SPACETREE_TIMEVARYING_H_

#include <list>

#include "ITL_util.h"
#include "ITL_histogram.h"
#include "ITL_partition_quadtree.h"
#include "ITL_partition_octree.h"
#include "ITL_distcomputer.h"

#include "FEL_util.h"
#include "FEL_range.h"
#include "FEL_domain.h"
#include "FEL_transformation.h"
#include "FEL_spacetreenode.h"

template <class T>
class FEL_spacetree_timevarying
{
	enum partitionType { QUADTREE = 1,
						 OCTREE = 2 };

	enum transformTypes { ROTATION = 0,
						  REFLECTION = 1 };

	// The list of all domains, needed for
	// comparison purpose while constructing the
	// range block tree.
	list< FEL_domain<T> >* domainList;
	int nDomainTotal;

	int globalRangeIndex;
	int nRangeTotal;
	int nLeafRange;

	float matchThreshold;
	float entropyThreshold;
	float blockSizeThreshold[3];

	FEL_transformation* trManager;

	ITL_partition_quadtree<T> quadtreePartitioner;
	ITL_partition_octree<T> octreePartitioner;

	ITL_histogram *histogram;
	int nBin;
	int distMeasureType;

	int nDim;
	int nTimeStep;
	//ITL_field_regular<T>* inputDataField;
	ITL_field_regular<int>** inputBinField;
	//ITL_field_regular<VECTOR3>* downsampledDataField;
	//ITL_field_regular<int>* downsampledDataBinField;

	int nChild;
	int partitionType;
	FEL_spacetreenode<T>* root;
	FEL_spacetreenode<T>* curNode;

	double* rangeLimitList;
	double* errorList;
	double* errorColorList;
	double* matchedDomainIdList;
	double* matchFreqList;

public:

	/**
	 * Constructor
	 */
	FEL_spacetree_timevarying()
	{
		nDim = 3;
		globalRangeIndex = 0;
		nRangeTotal = 0;
		nLeafRange = 0;
		matchThreshold = 0.01;
		//blockSizeThreshold[0] = blockSizeThreshold[1] = 16;
		distMeasureType = 0;

		root = NULL;
		curNode = NULL;
		histogram = NULL;
		inputBinField = NULL;
		//downsampledDataBinField = NULL;
		trManager = NULL;
	}

	/**
	 * Constructor
	 */
	void
	initialize( ITL_field_regular<int>** inputbinfield,
				list< FEL_domain<T> >* dList,
				ITL_histogram** hg, int nbin,
				int* rangeBlockSize,
				int distmeasuretype,
				int partitiontype,
				float matchthreshold,
				float entropythreshold,
			    FEL_transformation* tr )
	{
		inputBinField = inputbinfield;

		// Get domain list
		domainList = dList;
		nDomainTotal = domainList->size();

		nDim = 3;

		// Get histogram
		histogram = *hg;
		nBin = nbin;

		// Get transformation
		trManager = tr;

		// Set encoding parameters
		blockSizeThreshold[0] = rangeBlockSize[0];
		blockSizeThreshold[1] = rangeBlockSize[1];
		blockSizeThreshold[2] = rangeBlockSize[2];
		matchThreshold = matchthreshold;
		entropyThreshold = entropythreshold;
		partitionType = partitiontype;
		distMeasureType = distmeasuretype;

		// Initialize other tree properties
		if( partitionType == QUADTREE )	nChild = 4;
		else if( partitionType == OCTREE )	nChild = 8;
		globalRangeIndex = 0;
		nRangeTotal = 0;
		nLeafRange = 0;
		root = NULL;
		curNode = NULL;
		rangeLimitList = NULL;

		// Initialize matlab variables to NULL
		errorList = NULL;
		errorColorList = NULL;
		matchedDomainIdList = NULL;
		matchFreqList = NULL;

	}// End function

	/**
	 *
	 */
	void
	setMatlabVariables( double *rangelimits, double* errors,
						double* errorcolors, double* matcheddomains,
						double* matchedfrequencies )
	{
		rangeLimitList = rangelimits;
		errorList = errors;
		errorColorList = errorcolors;
		matchedDomainIdList = matcheddomains;
		matchFreqList = matchedfrequencies;

	}// End function

	/**
	 * Initialization of tree.
	 * Create a single node tree that contains
	 * the entire spatial domain of the field.
	 */
	void
	initTree( float* low, float* high )
	{
		// Create root node using the entire spatial domain
		root = new FEL_spacetreenode<T>( low, high );

		// Set properties
		root->setLevel( 0 );

		// Find best matching domain for this range block
		double matchError;
		fprintf( stderr, "Root node: <%g, %g, %g> <%g, %g, %g>\n", low[0], low[1], low[2], high[0], high[1], high[2] );

		// Set current node
		curNode = root;

	}// End function

	/**
	 * Recursive function that starts from the root and
	 * subdivides the spatial domain into range blocks.
	 */
	void
	createTree( FEL_spacetreenode<T>* cur,
				float* parentLow, float* parentHigh,
				int level )
	{
		int nDim = 3;

		double curRangeFreqList[nBin];
		int rotAmount = 0;
		bool tFlag = false;

		// Set properties of current node
		cur->setLevel( level );
		cur->setSpatialRange( parentLow,  parentHigh );

		fprintf( stderr, "Current range-block: <%g %g %g> <%g %g %g>\n",
							   parentLow[0], parentLow[1], parentLow[2],
							   parentHigh[0], parentHigh[1], parentHigh[2] );

		// Increase range block count
		nRangeTotal++;

		float childSize[4];
		float* childLow[8];
		float* childHigh[8];
		for( int i=0; i<8; i++ )
		{
			childLow[i] = new float[nDim];
			childHigh[i] = new float[nDim];
		}

		// If the range block size has fallen below some threshold, say <8,8,2>
		// // then stop recursion
		if( ( partitionType == QUADTREE &&
				((parentHigh[0]-parentLow[0]+1) <= blockSizeThreshold[0] ||
				 (parentHigh[1]-parentLow[1]+1) <= blockSizeThreshold[1]
				)
			)
			||
			( partitionType == OCTREE &&
				( (parentHigh[0]-parentLow[0]+1) <= blockSizeThreshold[0] ||
				  (parentHigh[1]-parentLow[1]+1) <= blockSizeThreshold[1] ||
				  (parentHigh[2]-parentLow[2]+1) <= blockSizeThreshold[2]
				)
			)
		  )
		  {
			// Store leaf rangeblock extent and other information for matlab plotting
			//rangeLimitList[globalRangeIndex*7] = parentLow[0];
			//rangeLimitList[globalRangeIndex*7+1] = parentHigh[0];
			//rangeLimitList[globalRangeIndex*7+2] = parentLow[1];
			//rangeLimitList[globalRangeIndex*7+3] = parentHigh[1];
			//rangeLimitList[globalRangeIndex*7+4] = parentLow[2];
			//rangeLimitList[globalRangeIndex*7+5] = parentHigh[2];
			//rangeLimitList[globalRangeIndex*7+6] = level;
			//errorList[globalRangeIndex] = matchError;
			//errorColorList[globalRangeIndex*3+2] = 0.3 + 0.7*(matchError / 100);
			//errorColorList[globalRangeIndex*3] = errorColorList[globalRangeIndex*3+1] = 0.005;
			//fprintf( stderr, "Stop recursing ...\n" );

			// Increase leaf range count
			globalRangeIndex++;
			nLeafRange ++;

			return;
		}
		// If the block is large enough and still no good match has been found
		// Subdivide and recurse
		else
		{
			if( partitionType == QUADTREE )
			{
				// Create partition of the current range block
				quadtreePartitioner.partition_Quadtree( parentLow, parentHigh, childLow, childHigh );

				// Create four children
				cur->setChildren( 4 );

				// Bottom left
				for( int i=0; i<3; i++ )
					childSize[i] = (int)(childHigh[0][i] - childLow[0][i] + 1);
				createTree( cur->getChild(0), childLow[0], childHigh[0], level+1 );

				// Bottom right
				for( int i=0; i<3; i++ )
					childSize[i] = (int)(childHigh[1][i] - childLow[1][i] + 1);
				createTree( cur->getChild(1), childLow[1], childHigh[1], level+1 );

				// Top right
				for( int i=0; i<3; i++ )
					childSize[i] = (int)(childHigh[2][i] - childLow[2][i] + 1);
				createTree( cur->getChild(2), childLow[2], childHigh[2], level+1 );

				// Top left
				for( int i=0; i<3; i++ )
					childSize[i] = (int)(childHigh[3][i] - childLow[3][i] + 1);
				createTree( cur->getChild(3), childLow[3], childHigh[3], level+1 );

			}// end if: partition type is quadtree
			else if( partitionType == OCTREE )
			{
				// Create partition of the current range block
				octreePartitioner.partition_Octree( parentLow, parentHigh, childLow, childHigh );

				// Create four children
				cur->setChildren( 8 );

				// 0
				for( int i=0; i<3; i++ )
					childSize[i] = (int)(childHigh[0][i] - childLow[0][i] + 1);
				createTree( cur->getChild(0), childLow[0], childHigh[0], level+1 );

				// 1
				for( int i=0; i<3; i++ )
					childSize[i] = (int)(childHigh[1][i] - childLow[1][i] + 1);
				createTree( cur->getChild(1), childLow[1], childHigh[1], level+1 );

				// 2
				for( int i=0; i<3; i++ )
					childSize[i] = (int)(childHigh[2][i] - childLow[2][i] + 1);
				createTree( cur->getChild(2), childLow[2], childHigh[2], level+1 );

				// 3
				for( int i=0; i<3; i++ )
					childSize[i] = (int)(childHigh[3][i] - childLow[3][i] + 1);
				createTree( cur->getChild(3), childLow[3], childHigh[3], level+1 );

				// 4
				for( int i=0; i<3; i++ )
					childSize[i] = (int)(childHigh[4][i] - childLow[4][i] + 1);
				createTree( cur->getChild(4), childLow[4], childHigh[4], level+1 );

				// 5
				for( int i=0; i<3; i++ )
					childSize[i] = (int)(childHigh[5][i] - childLow[5][i] + 1);
				createTree( cur->getChild(5), childLow[5], childHigh[5], level+1 );

				// 6
				for( int i=0; i<3; i++ )
					childSize[i] = (int)(childHigh[6][i] - childLow[6][i] + 1);
				createTree( cur->getChild(6), childLow[6], childHigh[6], level+1 );

				// 7
				for( int i=0; i<3; i++ )
					childSize[i] = (int)(childHigh[7][i] - childLow[7][i] + 1);
				createTree( cur->getChild(7), childLow[7], childHigh[7], level+1 );

			}// end if: partition type is octree

		}

		// Free resources
		for( int i=0; i<8; i++ )
		{
			delete [] childLow[i];
			delete [] childHigh[i];
		}


	}// End function

	/**
	 * Recursive function that starts from the root and
	 * for each tree node, find the matching domain
	 */
	void
	fillTree( FEL_spacetreenode<T> *cur, int t = 0 )
	{
		int nDim = 3;
		float curLow[3], curHigh[3];
		double matchError = 100;
		int nPointRange;
		int* rangeNodeData = NULL;
		double curRangeFreqList[nBin];
		double rangeTvFreqList[nBin*nTimeStep];
		//ITL_field_regular<int>** binFieldPtr = (*inputBinField);

		// Get properties of current node
		cur->getSpatialRange( curLow, curHigh );
		nPointRange = cur->getRangeBlock()->getSize();


		fprintf( stderr, "Current range-block: <%g %g %g> <%g %g %g>\n",
						 curLow[0], curLow[1], curLow[2],
						 curHigh[0], curHigh[1], curHigh[2] );

		// Get data
		rangeNodeData = new int[nPointRange];
		(*inputBinField)->getDataBetween( curLow, curHigh, rangeNodeData );

		// Get distribution
		cur->getRangeBlockDisrtibution( rangeTvFreqList, nBin, nTimeStep );

		// Convert to distribution
		FEL_util<T>::computeHistogramFromBinField( rangeNodeData, nPointRange, nBin, curRangeFreqList );

		// Update distribution for particular time step
		memcpy( rangeTvFreqList + t*nBin, curRangeFreqList, sizeof(double)*nBin );

		// Update distribution of range block in the tree
		cur->setRangeBlockDisrtibution( rangeTvFreqList, nBin, nTimeStep );

		if( cur->isLeaf() )
			return;
		else
		{
			//fprintf( stderr, "recursing down: <%g,%g>\n", lim[0], lim[1] );
			for( int iC=0; iC<nChild; iC++ )
				fillTree( cur->getChild(iC), t );

		}// end inner else

	}// End function

	/**
	 * Recursive function that starts from the root and
	 * for each tree node, find the matching domain
	 */
	void
	matchTree( FEL_spacetreenode<T> *cur )
	{
		int nDim = 3;
		float curLow[3], curHigh[3];
		float curEntropy;
		double matchError = 100;
		int rotAmount;
		bool tFlag;
		int nPointRange;
		int* rangeNodeData = NULL;
		double curRangeTvFreqList[nBin*nTimeStep];

		// Get limits of current node
		cur->getSpatialRange( curLow,  curHigh );

		fprintf( stderr, "Current range-block: <%g %g %g> <%g %g %g>\n",
						 curLow[0], curLow[1], curLow[2],
						 curHigh[0], curHigh[1], curHigh[2] );

		// Get distribution of the current block
		cur->getRangeBlockDisrtibution( curRangeTvFreqList, nBin, nTimeStep );

		// Find best matching domain for this range block
		int bestMatchedDomainId = matchRangeToDomains2( curRangeTvFreqList,
														&matchError,
														&rotAmount, &tFlag );

		curEntropy = FEL_util<double>::computeEntropy_TimeHistogram( curRangeTvFreqList, nBin, nTimeStep, false );

		fprintf( stderr, "Entropy: %g Match error: %g Match Transform: %d %d Match ID: %d\n",
						 curEntropy, matchError, rotAmount, (int)tFlag, bestMatchedDomainId );

		//cur->setEntropy( curEntropy );

		// Add range block to the tree node
		// Store range distribution
		// Store matching domain information as well
		cur->setRangeBlock( nRangeTotal,
							curLow, curHigh,
							curRangeTvFreqList,
							nBin, nTimeStep,
							curEntropy,
							rotAmount, tFlag,
							matchError, bestMatchedDomainId );


		// Add matching range information to the domain
		typename list<FEL_domain<T> >::iterator tmpIter = domainList->begin();
		for( int iD=0; iD<bestMatchedDomainId; iD++ )
			++tmpIter;
		tmpIter->addMatchingRange( cur->getRangeBlock()->getID(), matchError, rotAmount, tFlag  );

		if( cur->isLeaf() )
			return;
		else
		{
			//fprintf( stderr, "recursing down: <%g,%g>\n", lim[0], lim[1] );
			for( int iC=0; iC<nChild; iC++ )
				matchTree( cur->getChild(iC) );

		}// end inner else

	}// End function

	list< FEL_range<T> >
	searchTreeByDomainID( FEL_spacetreenode<T> *cur, int qDomainID )
	{
		list< FEL_range<T> > localSearchResultList;

		// Get domain ID associated to this node
		int curNodeMatchingDomain = cur->getRangeBlock()->getMatchedDomainId();
		//printf( "At Node: <%g,%g>\n", lim[0], lim[1] );

		// If the node's matching domain ID matches the query
		// domain ID, then return current range
		if( qDomainID == curNodeMatchingDomain )
		{
			float bounds[6];
			cur->getSpatialRange( bounds );
			//printf( "Returned range: <%g,%g,%g>, <%g,%g,%g>\n", bounds[0], bounds[1], bounds[2],
			//		 	 	 	 	 	 	 	 	 	 	 	bounds[3], bounds[4], bounds[5] );

			localSearchResultList.push_back( *(cur->getRangeBlock()) );
			return localSearchResultList;
		}
		else
		{
			if( cur->isLeaf() )
				return localSearchResultList;
			else
			{
				//printf( "recursing down: <%g,%g>\n", lim[0], lim[1] );
				list< FEL_range<T> > cList;
				for( int iC=0; iC<nChild; iC++ )
				{
					cList.clear();
					cList = searchTreeByDomainID( cur->getChild(iC), qDomainID );
					typename list< FEL_range<T> >::iterator it = localSearchResultList.end();
					if( !cList.empty() )
						localSearchResultList.insert( it, cList.begin(), cList.end() );
				}// end for
			}// end inner else
		}// end outer else

		return localSearchResultList;

	}// End function

	list< FEL_range<T>* >
	searchTreeByLevel( FEL_spacetreenode<T> *cur, int qLevel )
	{
		list< FEL_range<T>* > localSearchResultList;

		// Get domain ID associated to this node
		int curNodeLevel = cur->getLevel();
		//printf( "At Node: <%g,%g>\n", lim[0], lim[1] );

		// If the node's level matches the query
		// level, then return the nodes corresponding range block
		if( qLevel == curNodeLevel )
		{
			//float bounds[6];
			//cur->getSpatialRange( bounds );
			//printf( "Returned range: <%g,%g,%g>, <%g,%g,%g>\n", bounds[0], bounds[1], bounds[2],
			//		 	 	 	 	 	 	 	 	 	 	 	bounds[3], bounds[4], bounds[5] );

			localSearchResultList.push_back( cur->getRangeBlock() );
			return localSearchResultList;
		}
		else
		{
			if( cur->isLeaf() )
				return localSearchResultList;
			else
			{
				//printf( "recursing down: <%g,%g>\n", lim[0], lim[1] );
				list< FEL_range<T>* > cList;
				for( int iC=0; iC<nChild; iC++ )
				{
					cList.clear();
					cList = searchTreeByLevel( cur->getChild(iC), qLevel );
					typename list< FEL_range<T>* >::iterator it = localSearchResultList.end();
					if( !cList.empty() )
						localSearchResultList.insert( it, cList.begin(), cList.end() );
				}// end for
			}// end inner else
		}// end outer else

		return localSearchResultList;

	}// End function

	FEL_range<T>*
	searchTreeByRangeID( FEL_spacetreenode<T> *cur, int qRangeID )
	{
		// Get range ID associated to this node
		int curNodeID = cur->getRangeBlock()->getID();
		//printf( "At Node: <%d,%d>\n", curNodeID, cur->getLevel() );

		// If the node's ID matches the query
		// range ID, then return current range
		if( qRangeID == curNodeID )
		{
			float bounds[6];
			cur->getSpatialRange( bounds );
			//printf( "Returned range: <%g,%g,%g>, <%g,%g,%g>\n", bounds[0], bounds[1], bounds[2],
			//		 	 	 	 	 	 	 	 	 	 	 	bounds[3], bounds[4], bounds[5] );
			assert( cur->getRangeBlock() != NULL );
			return cur->getRangeBlock();
		}
		else
		{
			if( cur->isLeaf() )
				return NULL;
			else
			{
				int iC=0;
				FEL_range<T>* ret = NULL;
				while( iC<nChild )
				{
					ret = searchTreeByRangeID( cur->getChild(iC), qRangeID );

					if( ret != NULL )
						return ret;

					iC++;

				}// end while
			}// end inner else
		}// end outer else

		return NULL;

	}// End function

	list< FEL_range<T>* >
	getLeafNodes( FEL_spacetreenode<T> *cur )
	{
		list< FEL_range<T>* > localSearchResultList;

		if( cur->isLeaf() )
		{
			localSearchResultList.push_back( cur->getRangeBlock() );
			return localSearchResultList;
		}
		else
		{
			list< FEL_range<T>* > cList;
			for( int iC=0; iC<nChild; iC++ )
			{
				cList.clear();
				cList = getLeafNodes( cur->getChild(iC) );
				typename list< FEL_range<T>* >::iterator it = localSearchResultList.end();
				if( !cList.empty() )
					localSearchResultList.insert( it, cList.begin(), cList.end() );

			}// end for

		}// end outer else

		return localSearchResultList;
	}

	list< FEL_range<T> >
	searchLeafNodesByRegion( FEL_spacetreenode<T>* cur, float* queryBounds )
	{
		float curRangeBounds[6];
		list< FEL_range<T> > localSearchResultList;
		list< FEL_range<T> > cList;
		typename list< FEL_range<T> >::iterator iter;
		//list< FEL_range<T>* > localSearchResultList;
		//list< FEL_range<T>* > cList;
		///typename list< FEL_range<T>* >::iterator iter;
		float queryBoundsLocal[6];

		// Get bounds for this node
		cur->getSpatialRange( curRangeBounds );

		#ifdef DEBUG_MODE
		fprintf( stderr, "Current Node: %d <%g %g %g>-<%g %g %g>\n",
				 cur->getLevel(),
				 curRangeBounds[0], curRangeBounds[1], curRangeBounds[2],
				 curRangeBounds[3], curRangeBounds[4], curRangeBounds[5]);
		#endif

		// Check for overlap
		// If the current node does not overlap,
		// none of its children will, so recursion can be stopped.
		// If the current node does overlap,
		// and it is a leaf node, then return, else continue redursion.
		if( !FEL_util<float>::isOverlapping( queryBounds, curRangeBounds ) )
			return localSearchResultList;
		else
		{
			if( cur->isLeaf() )
			{
				// Add leaf node and Return result
				#ifdef DEBUG_MODE
				fprintf( stderr, "At leaf node ..\n" );
				#endif
				localSearchResultList.push_back( *cur->getRangeBlock() );

				return localSearchResultList;
			}
			else
			{
				//fprintf( stderr, "Non-leaf node, hence recursing down ...\n" );
				//memcpy( queryBoundsLocal, (float*)queryBounds, sizeof(float)*6 );
				for( int i=0; i<6; i++ )
					queryBoundsLocal[i] = queryBounds[i];

				for( int iC=0; iC<nChild; iC++ )
				{
					cList.clear();
					cList = searchLeafNodesByRegion( cur->getChild(iC), queryBoundsLocal );
					if( cList.size() == 0 )
						continue;
					iter = localSearchResultList.end();
					if( !cList.empty() )
						localSearchResultList.insert( iter, cList.begin(), cList.end() );
				}// end for

			} // end inner else : leaf node or not

		}// end outer else : has overlap or not

		return localSearchResultList;

	}

	void
	setMinimumBlockSize( float* blocksize )
	{
		for( int i=0; i<3; i++ )
			blockSizeThreshold[i] = blocksize[i];
	}// End function

	int
	matchRangeToDomains2( double* rangeFreqList,
						  double *error,
						  int* rotateAmount,
						  bool* isReflected
						 )
	{
		int bestMatchedDomainId = -1;
		double domainFreqList[nBin];
		double e = 100;
		double minError = 100;
		int rotAmount = 0;
		bool isRef = false;

		// Scan through all domains to find matching domain
		int domainID = 0;
		for( typename list< FEL_domain<T> >::iterator domainIter = domainList->begin();
			 domainIter != domainList->end();
			 domainIter++ )
		{
			int bounds[6];

			domainIter->getDistribution( domainFreqList, &nBin );
			#ifdef DEBUG_MODE
			double m = ITL_util<double>::Min( domainFreqList, nBin );
			double M = ITL_util<double>::Max( domainFreqList, nBin );
			fprintf( stderr, "%d-th Domain distribution min/max: %g, %g\n", domainID, m, M );
			#endif

			// Transform domains by data
			//e =  matchRangeToDomain2( rangeFreqList, *domainIter, optTransformationMatrix );
			// Transform domains by distribution
			e =  FEL_util<double>::matchRangeToDomainByDistribution( distMeasureType,
																	 rangeFreqList, domainFreqList,
																	 nBin,
																	 &rotAmount, &isRef );

			if( e < minError )
			{
				//cout << "-> " << rotAmount << " " << isRef << endl;
				minError = e;
		  	    (*rotateAmount) = rotAmount;
				(*isReflected) = isRef;
				bestMatchedDomainId = domainID;
			}

			domainID ++;
		}

		*error = minError;
		return bestMatchedDomainId;

	}// End function

	void
	printDomains()
	{
		int domainBounds[6];

		int domainID = 0;
		for( typename list< FEL_domain<T> >::iterator domainIter = domainList->begin();
			 domainIter != domainList->end(); domainIter++ )
		{
			domainIter->getBounds( domainBounds );
			printf( "%d-th domain: <%d %d %d> <%d %d %d>, Entropy: %g\n", domainID,
									domainBounds[0], domainBounds[1], domainBounds[2],
									domainBounds[3], domainBounds[4], domainBounds[5],
									domainIter->getEntropy() );
			domainID ++;
		}// End for

	}// End function

	int
	countNumberOfTreeNodes( FEL_spacetreenode<T> *cur )
	{
		if( cur->isLeaf() )
			return 1;

		int count = 0;
		for( int iC=0; iC<nChild; iC++ )
			count += countNumberOfTreeNodes( cur->getChild(iC) );

		return (count+1);
	}// End function

	int
	countNumberOfLeafNodes( FEL_spacetreenode<T> *cur )
	{
		if( cur->isLeaf() )
			return 1;

		int count = 0;
		for( int iC=0; iC<nChild; iC++ )
			count += countNumberOfTreeNodes( cur->getChild(iC) );

		return count;

	}// End function

	FEL_spacetreenode<T>*
	getRoot()
	{
		return root;
	}

	int
	getNumRangeBlock()
	{
		return nRangeTotal;
	}

	int
	getNumLeafRangeBlock()
	{
		return nLeafRange;
	}

	void
	setBinField( ITL_field_regular<int>* fieldptr )
	{
		// Do not delete because it is a shallow copy
		inputBinField = fieldptr;
	}

	void
	printTreeNodes()
	{

	}// End function

	void
	printLeafNodes( FEL_spacetreenode<T> *cur )
	{
		//if( cur->isLeaf() )
		//{
		//	cur->pr
		//}

		//int count = 0;
		//for( int iC=0; iC<nChild; iC++ )
		//	count += countNumberOfTreeNodes( cur->getChild(iC) );


	}// End function
};

#endif
/* FEL_SPACETREE_TIMEVARYING_H_ */
