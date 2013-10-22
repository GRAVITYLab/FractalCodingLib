/*
 * FEL_spacetree.h
 *
 *  Created on: Feb 23, 2012
 *      Author: abon
 */

#ifndef FEL_SPACETREE_H_
#define FEL_SPACETREE_H_

#include <list>

#include "ITL_util.h"
#include "ITL_histogram.h"
#include "ITL_partition_quadtree.h"
#include "ITL_partition_octree.h"
#include "ITL_distcomputer.h"

#include "FEL_util.h"
#include "FEL_range.h"
#include "FEL_domain.h"
#include "FEL_matlab.h"
#include "FEL_transformation.h"
#include "FEL_spacetreenode.h"

#define EPSILON 0.01

template <class T>
class FEL_spacetree
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
	//ITL_field_regular<T>* inputDataField;
	ITL_field_regular<int>* inputBinField;
	//ITL_field_regular<VECTOR3>* downsampledDataField;
	//ITL_field_regular<int>* downsampledDataBinField;

	int nChild;
	int partitionType;
	FEL_spacetreenode<T>* root;
	FEL_spacetreenode<T>* curNode;

	// Plotting
	bool matlab_plotting_on;
	FEL_matlab** plotter;

public:

	int noGoodMatchFound;

	/**
	 * Constructor
	 */
	FEL_spacetree()
	{
		nDim = 3;
		globalRangeIndex = 0;
		nRangeTotal = 0;
		nLeafRange = 0;
		matchThreshold = 0.01;
		//blockSizeThreshold[0] = blockSizeThreshold[1] = 16;
		distMeasureType = 0;
		matlab_plotting_on = false;

		root = NULL;
		curNode = NULL;
		histogram = NULL;
		inputBinField = NULL;
		//downsampledDataBinField = NULL;
		trManager = NULL;
		noGoodMatchFound = 0;

	}

	/**
	 * Constructor
	 */
	void
	initialize( list< FEL_domain<T> >* dList,
				ITL_histogram** hg, int nbin,
				//ITL_field_regular<T>** dataFieldPtr,
			    ITL_field_regular<int>** binFieldPtr,
				int* rangeBlockSize,
				int distmeasuretype,
				int partitiontype,
				float matchthreshold,
				float entropythreshold,
			    FEL_transformation* tr )
	{
		// Get domain list
		domainList = dList;
		nDomainTotal = domainList->size();

		// Get Data
		//inputDataField = *dataFieldPtr;
		inputBinField = *binFieldPtr;
	//	downsampledDataField = *dsampledDataFieldPtr;
	//	downsampledDataBinField = *dsampledDataBinFieldPtr;
		nDim = inputBinField->getNumDim();

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
		matlab_plotting_on = false;

		// Initialize other tree properties
		if( partitionType == QUADTREE )	nChild = 4;
		else if( partitionType == OCTREE )	nChild = 8;
		globalRangeIndex = 0;
		nRangeTotal = 0;
		nLeafRange = 0;
		root = NULL;
		curNode = NULL;
		noGoodMatchFound = 0;

	}// End function

	/**
	 *
	 */
	void
	setMatlabPlotter( FEL_matlab** p, bool flag )
	{
		plotter = p;
		matlab_plotting_on = flag;
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
		MATRIX3 optimalTrMatrix( VECTOR3( 1, 0, 0 ),
								 VECTOR3( 0, 1, 0 ),
								 VECTOR3( 0, 0, 1 ) );
		//int bestMatchedDomainId = matchRangeToDomains2( low, high, &matchError, &optimalTrMatrix );
		fprintf( stderr, "Root node: <%g, %g, %g> <%g, %g, %g>\n", low[0], low[1], low[2], high[0], high[1], high[2] );
		//fprintf( stderr, "Root node: Match error: %g Match ID: %d\n", matchError, bestMatchedDomainId );

		// Set current node
		curNode = root;

	}// End function

	/**
	 * Recursive function that starts from the root and
	 * subdivides the spatial domain into range blocks.
	 */
	void
	createTree2( FEL_spacetreenode<T> *cur,
				float* parentLow, float* parentHigh,
				int level )
	{
		int nDim = 3;
		double matchError = 100;
		double muError = 100;
		double sdError = 100;
		double curRangeFreqList[nBin];
		float curEntropy;
		int rotAmount = 0;
		bool tFlag = false;

		bool stopPartition = true;

		// Set properties of current node
		cur->setLevel( level );
		cur->setSpatialRange( parentLow,  parentHigh );

		#ifdef DEBUG_MODE
		fprintf( stderr, "Current range-block: <%g %g %g> <%g %g %g>\n",
							   parentLow[0], parentLow[1], parentLow[2],
							   parentHigh[0], parentHigh[1], parentHigh[2] );
		#endif

		// Find best matching domain for this range block
		int bestMatchedDomainId = matchRangeToDomains2( parentLow, parentHigh,
														&matchError,
														//&optimalTrMatrix,
														&rotAmount, &tFlag,
														curRangeFreqList,
														&muError, &sdError );
		curEntropy = ITL_entropycore::computeEntropy_HistogramBased2( curRangeFreqList, nBin, false );

		#ifdef DEBUG_MODE
		fprintf( stderr, "Entropy: %g Match error: %g Match Transform: %d %d Match ID: %d\n",
						 curEntropy, matchError, rotAmount, (int)tFlag, bestMatchedDomainId );
		#endif

		//cur->setEntropy( curEntropy );

		// Add range block to the tree node
		// Store range distribution
		// Store matching domain information as well
		cur->setRangeBlock( nRangeTotal,
							parentLow, parentHigh,
							curRangeFreqList, nBin, 1,
							curEntropy,
							rotAmount, tFlag,
							matchError, bestMatchedDomainId );

		// Add matching range information to the domain
		typename list<FEL_domain<T> >::iterator tmpIter = domainList->begin();
		for( int iD=0; iD<bestMatchedDomainId; iD++ )
			++tmpIter;
		tmpIter->addMatchingRange( cur->getRangeBlock()->getID(), matchError, rotAmount, tFlag  );

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
		// OR the matching error has fallen below the threshold
		// then stop recursion
		if( matchError > matchThreshold ) stopPartition = false;

		if( curEntropy > entropyThreshold ) stopPartition = false;

		if( ( partitionType == QUADTREE &&
			  ( (parentHigh[0]-parentLow[0]+1) <= blockSizeThreshold[0] ||
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
			stopPartition = true;

		if( stopPartition == true )
		{
			if( matchError > matchThreshold )
			{
				noGoodMatchFound ++;

				// Mark the range (the range histogram is being stored anyway)
				cur->setRangeMatchQualityFlag( false );


			}


			// Store leaf rangeblock extent and other information for matlab plotting
			if( matlab_plotting_on )
			{
				(*plotter)->rangeLimitList[globalRangeIndex*7] = parentLow[0];
				(*plotter)->rangeLimitList[globalRangeIndex*7+1] = parentHigh[0];
				(*plotter)->rangeLimitList[globalRangeIndex*7+2] = parentLow[1];
				(*plotter)->rangeLimitList[globalRangeIndex*7+3] = parentHigh[1];
				(*plotter)->rangeLimitList[globalRangeIndex*7+4] = parentLow[2];
				(*plotter)->rangeLimitList[globalRangeIndex*7+5] = parentHigh[2];
				(*plotter)->rangeLimitList[globalRangeIndex*7+6] = level;

				if( matchError > matchThreshold )
				{
					(*plotter)->unMatchedRangeFlagList[globalRangeIndex] = 1;
					(*plotter)->errorList[globalRangeIndex] = 0;
					(*plotter)->muErrorList[globalRangeIndex] = 0;
					(*plotter)->sdErrorList[globalRangeIndex] = 0;
				}
				else
				{
					(*plotter)->errorList[globalRangeIndex] = matchError;
					(*plotter)->muErrorList[globalRangeIndex] = muError;
					(*plotter)->sdErrorList[globalRangeIndex] = sdError;
				}
				(*plotter)->rangeEntropyList[globalRangeIndex] = curEntropy;
				(*plotter)->matchedDomainIdList[globalRangeIndex] = bestMatchedDomainId;
				if( bestMatchedDomainId != -1 )
					(*plotter)->domainMatchFreqList[bestMatchedDomainId] ++;
			}

			#ifdef DEBUG_MODE
			fprintf( stderr, "%d-th range maps to %d-th domain with %g error ...\n", globalRangeIndex,
																	bestMatchedDomainId, matchError );
			#endif
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
				createTree2( cur->getChild(0), childLow[0], childHigh[0], level+1 );

				// Bottom right
				for( int i=0; i<3; i++ )
					childSize[i] = (int)(childHigh[1][i] - childLow[1][i] + 1);
				createTree2( cur->getChild(1), childLow[1], childHigh[1], level+1 );

				// Top right
				for( int i=0; i<3; i++ )
					childSize[i] = (int)(childHigh[2][i] - childLow[2][i] + 1);
				createTree2( cur->getChild(2), childLow[2], childHigh[2], level+1 );

				// Top left
				for( int i=0; i<3; i++ )
					childSize[i] = (int)(childHigh[3][i] - childLow[3][i] + 1);
				createTree2( cur->getChild(3), childLow[3], childHigh[3], level+1 );

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
				createTree2( cur->getChild(0), childLow[0], childHigh[0], level+1 );

				// 1
				for( int i=0; i<3; i++ )
					childSize[i] = (int)(childHigh[1][i] - childLow[1][i] + 1);
				createTree2( cur->getChild(1), childLow[1], childHigh[1], level+1 );

				// 2
				for( int i=0; i<3; i++ )
					childSize[i] = (int)(childHigh[2][i] - childLow[2][i] + 1);
				createTree2( cur->getChild(2), childLow[2], childHigh[2], level+1 );

				// 3
				for( int i=0; i<3; i++ )
					childSize[i] = (int)(childHigh[3][i] - childLow[3][i] + 1);
				createTree2( cur->getChild(3), childLow[3], childHigh[3], level+1 );

				// 4
				for( int i=0; i<3; i++ )
					childSize[i] = (int)(childHigh[4][i] - childLow[4][i] + 1);
				createTree2( cur->getChild(4), childLow[4], childHigh[4], level+1 );

				// 5
				for( int i=0; i<3; i++ )
					childSize[i] = (int)(childHigh[5][i] - childLow[5][i] + 1);
				createTree2( cur->getChild(5), childLow[5], childHigh[5], level+1 );

				// 6
				for( int i=0; i<3; i++ )
					childSize[i] = (int)(childHigh[6][i] - childLow[6][i] + 1);
				createTree2( cur->getChild(6), childLow[6], childHigh[6], level+1 );

				// 7
				for( int i=0; i<3; i++ )
					childSize[i] = (int)(childHigh[7][i] - childLow[7][i] + 1);
				createTree2( cur->getChild(7), childLow[7], childHigh[7], level+1 );

			}// end if: partition type is octree

		}

		// Free resources
		for( int i=0; i<8; i++ )
		{
			delete [] childLow[i];
			delete [] childHigh[i];
		}


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
	matchRangeToDomains2( float *rangeLow, float *rangeHigh,
						  double *error,
						  //MATRIX3 *optTransformationMatrix,
						  int* rotateAmount,
						  bool* isReflected,
						  double* rangeFreqList,
						  double* muerror = NULL, double* sderror = NULL
						  )
	{
		int bestMatchedDomainId = -1;
		double domainFreqList[nBin];
		double e = 100;
		double minError = 100;
		int rotAmount = 0;
		bool isRef = false;
		double bestDomainFreqList[nBin];

		// Collect distribution for this range
		//
		// Get range block size
		int nPointRange = 1;
		for( int i=0; i<nDim; i++ )
			nPointRange *= (int)( rangeHigh[i] - rangeLow[i] + 1 );

		// Get range block data
		int* rangeBinData = new int[nPointRange];
		//int rangeBinData[nPointRange];
		inputBinField->getDataBetween( rangeLow, rangeHigh, rangeBinData );

		// Convert to distribution
		FEL_util<double>::computeHistogramFromBinField( rangeBinData, nPointRange, nBin, rangeFreqList );
		#ifdef DEBUG_MODE
		double mr = ITL_util<double>::Min( rangeFreqList, nBin );
		double Mr = ITL_util<double>::Max( rangeFreqList, nBin );
		fprintf( stderr, "Current range distribution min/max: %g, %g\n", mr, Mr );
		#endif

		// Compute range entropy
		//(*rangeEntropy) = (float)ITL_entropycore::computeEntropy_HistogramBased2( rangeFreqList, nBin, false );

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
				memcpy( bestDomainFreqList, domainFreqList, sizeof( double )*nBin );
			}

			domainID ++;
		}

		delete [] rangeBinData;

		*error = minError;

		// Optional analysis step: compute mean diffference and sd difference
		if( muerror != NULL )
			(*muerror) = FEL_util<double>::computeHistogramMeanError( rangeFreqList, bestDomainFreqList, nBin );
		if( sderror != NULL )
			(*sderror) = FEL_util<double>::computeHistogramSDError( rangeFreqList, bestDomainFreqList, nBin );

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

};

#endif
/* FEL_SPACETREE_H_ */
