/*
 * FEL_matlab.cpp
 *
 *  Created on: Mar 5, 2012
 *      Author: abon
 */

#include "FEL_matlab.h"

FEL_matlab::FEL_matlab()
{
	memset( domainMatchFreqList, 0, sizeof(int)*NMAX );
	memset( unMatchedRangeFlagList, 0, sizeof(int)*NMAX );

	engine = NULL;

	mxNBin = NULL;
	mxNTotalDomain = NULL;
	mxNTotalRange = NULL;

	mxNUnmatchedRange = NULL;
	mxUnmatchedRangeFlagList = NULL;

	mxErrorList = NULL;
	mxMuErrorList = NULL;
	mxSdErrorList = NULL;

	mxMatchedDomainList = NULL;
	mxDomainMatchFreqList = NULL;

	mxRangeLimitList = NULL;

	mxDomEntropyList = NULL;
	mxRangeEntropyList = NULL;

}// End constructor


FEL_matlab::~FEL_matlab()
{
	// Free memory
	mxDestroyArray( mxNBin );
	mxDestroyArray( mxNTotalRange );
	mxDestroyArray( mxNUnmatchedRange );
	mxDestroyArray( mxRangeLimitList );
	mxDestroyArray( mxUnmatchedRangeFlagList );
	mxDestroyArray( mxMatchedDomainList );
	mxDestroyArray( mxNTotalDomain );
	mxDestroyArray( mxErrorList );
	mxDestroyArray( mxMuErrorList );
	mxDestroyArray( mxSdErrorList );
	mxDestroyArray( mxDomEntropyList );

}// End function


void
FEL_matlab::initEngine()
{
	// Call engOpen with a NULL string. This starts a MATLAB process
    // on the current host using the command "matlab".
	fprintf( stderr, "Initializing Matlab engine ...\n" );

	if (!(engine = engOpen("")))
	{
		fprintf(stderr, "\nCan't start MATLAB engine\n");
		exit(1);
	}
}// End function

void
FEL_matlab::closeEngine()
{
	// Close matlab engine
	engClose( engine );

}// End function

void
FEL_matlab::createMatlabVariables()
{
	double tmp = nBin;
	mxNBin = mxCreateDoubleMatrix( 1, 1, mxREAL );
	memcpy( (void *)mxGetPr(mxNBin), (void *)(&tmp), sizeof(double) );

	tmp = nDomain;
	mxNTotalDomain = mxCreateDoubleMatrix( 1, 1, mxREAL );
	memcpy( (void *)mxGetPr(mxNTotalDomain), (void *)(&tmp), sizeof(double) );

	tmp = nRange;
	mxNTotalRange = mxCreateDoubleMatrix( 1, 1, mxREAL );
	memcpy( (void *)mxGetPr(mxNTotalRange), (void *)(&tmp), sizeof(double) );

	tmp = nUnmatchedRange;
	mxNUnmatchedRange = mxCreateDoubleMatrix( 1, 1, mxREAL );
	memcpy( (void *)mxGetPr( mxNUnmatchedRange ), (void *)(&tmp), sizeof(double) );

	mxUnmatchedRangeFlagList = mxCreateDoubleMatrix( 1, nRange, mxREAL );
	memcpy( (void *)mxGetPr(mxUnmatchedRangeFlagList), (void *)unMatchedRangeFlagList, sizeof(double)*(nRange) );

	mxDomainMatchFreqList = mxCreateDoubleMatrix( 1, nDomain, mxREAL );
	memcpy( (void *)mxGetPr( mxDomainMatchFreqList ), (void *)domainMatchFreqList, sizeof(double)*nDomain );

	mxRangeLimitList = mxCreateDoubleMatrix( 7, nRange, mxREAL );
	memcpy( (void *)mxGetPr(mxRangeLimitList), (void *)rangeLimitList, sizeof(double)*(nRange*7) );

	mxMatchedDomainList = mxCreateDoubleMatrix( 1, nRange, mxREAL );
	memcpy( (void *)mxGetPr(mxMatchedDomainList), (void *)matchedDomainIdList, sizeof(double)*nRange );

	mxErrorList = mxCreateDoubleMatrix( 1, nRange, mxREAL );
	memcpy( (void *)mxGetPr( mxErrorList ), (void *)errorList, sizeof(double)*nRange );

	mxMuErrorList = mxCreateDoubleMatrix( 1, nRange, mxREAL );
	memcpy( (void *)mxGetPr( mxMuErrorList ), (void *)muErrorList, sizeof(double)*nRange );

	mxSdErrorList = mxCreateDoubleMatrix( 1, nRange, mxREAL );
	memcpy( (void *)mxGetPr( mxSdErrorList ), (void *)sdErrorList, sizeof(double)*nRange );

	mxDomEntropyList = mxCreateDoubleMatrix( 1, nDomain, mxREAL );
	memcpy( (void *)mxGetPr( mxDomEntropyList ), (void *)domainEntropyList, sizeof(double)*nDomain );

	mxRangeEntropyList = mxCreateDoubleMatrix( 1, nRange, mxREAL );
	memcpy( (void *)mxGetPr( mxRangeEntropyList ), (void *)rangeEntropyList, sizeof(double)*nRange );

	assert( engine != NULL );
	// Place the variables into the MATLAB workspace
	engPutVariable( engine, "mxNBin", mxNBin );
	engPutVariable( engine, "mxNTotalDomain", mxNTotalDomain );
	engPutVariable( engine, "mxNTotalRange", mxNTotalRange );
	engPutVariable( engine, "mxNUnmatchedRange", mxNUnmatchedRange );
	engPutVariable( engine, "mxUnmatchedRangeFlagList", mxUnmatchedRangeFlagList );
	engPutVariable( engine, "mxDomainEntropyList", mxDomEntropyList );
	engPutVariable( engine, "mxRangeEntropyList", mxRangeEntropyList );
	engPutVariable( engine, "mxErrorList", mxErrorList );
	engPutVariable( engine, "mxMuErrorList", mxMuErrorList );
	engPutVariable( engine, "mxSdErrorList", mxSdErrorList );
	engPutVariable( engine, "mxRangeLimitList", mxRangeLimitList );
	engPutVariable( engine, "mxMatchedDomainList", mxMatchedDomainList );
	engPutVariable( engine, "mxDomainMatchFreqList", mxDomainMatchFreqList );

}

void
FEL_matlab::displayDomainUtilization_Matlab()
{
	fprintf( stderr, "Displaying domain utilization plot ...\n" );

	// Plot domain ID vs frequency
	engEvalString( engine, "displayDomainUtilization( mxDomainMatchFreqList );" );
	//engEvalString( engine, "figure" );
	//engEvalString( engine, "plot( mxDomainMatchFreqList, 'b', 'LineWidth', 2 );" );
	//engEvalString( engine, "xlabel('Domain ID');" );
	//engEvalString( engine, "ylabel( 'Match Frequency');" );
}

void
FEL_matlab::displayMatchingPairEntropyLevel_Matlab()
{
	fprintf( stderr, "Displaying entropy levels of range-domain pair ...\n" );

	// Plot range entropy vs matching domain entropy
	engEvalString( engine, "displayMatchingPairEntropyLevel( mxNTotalRange, mxRangeEntropyList, mxDomainEntropyList, mxMatchedDomainList );" );
	//engEvalString( engine, "figure" );
	//engEvalString( engine, "mxMatchedfDomainEntropyList = zeros( mxNRange, 1);" );
	//engEvalString( engine, "mxMatchedfDomainEntropyList = mxDomainEntropyList( mxMatchedDomainList );" );
	//engEvalString( engine, "L = [mxRangeEntropyList'; mxDomainEntropyList']';" );
	//engEvalString( engine, "L = sortrows(L);" );
	//engEvalString( engine, "plot( L(1,:), L(2:), 'r.', 'LineWidth', 2 );" );
	//engEvalString( engine, "xlabel('Range Entropy');" );
	//engEvalString( engine, "ylabel('Matching Domain Entropy');" );
}

void
FEL_matlab::displayMatchErrorDistribution_Matlab()
{
	fprintf( stderr, "Displaying error distribution ...\n" );

	// Plot the error as a histogram
	engEvalString( engine, "displayMatchErrorDistribution( mxErrorList, mxMuErrorList, mxSdErrorList );" );
	//engEvalString( engine, "figure;" );
	//engEvalString( engine, "hist( mxErrorList, 20 );" );
	//engEvalString( engine, "xlabel('Range-domain Match Error(%)');" );
	//engEvalString( engine, "ylabel('Range frequency');" );

	// Plot range size vs. error
	engEvalString( engine, "displayRangeEntropyVsError( mxRangeEntropyList, mxErrorList );" );
	//engEvalString( engine, "figure;" );
	//engEvalString( engine, "L = [mxRangeEntropyList'; mxErrorList']';" );
	//engEvalString( engine, "L = sortrows(L);" );
	//engEvalString( engine, "plot( L(1,:), L(2,:), 'b', 'LineWidth', 2 );" );
	//engEvalString( engine, "xlabel('Range Size');" );
	//engEvalString( engine, "ylabel('Error(%)');" );

}

void
FEL_matlab::displayUnMatchedRanges_Matlab()
{
	fprintf( stderr, "Displaying unmatched ranges ...\n" );
	// Draw the unmatched ranges
	engEvalString( engine, "displayUnMatchedRanges( mxRangeLimitList, mxNUnmatchedRange, mxUnmatchedRangeFlagList );" );

}

void
FEL_matlab::displayBlockwiseMatchErrorPlot_Matlab()
{

}// End function

void
FEL_matlab::displayDomainEntropies_Matlab()
{
	// Plot the error for each domain
	engEvalString( engine, "figure" );
	engEvalString( engine, "plot( sort( mxDomEntropyList ), 'r*' );" );
	engEvalString( engine, "title('Entropy of domains');" );
	engEvalString( engine, "xlabel('Domain ID');" );
	engEvalString( engine, "ylabel('Entropy');" );
	engEvalString( engine, "axis( [0 mxNTotalDomain-1 0 log2(mxNBin)] );" );

}// End function

/*
void
FEL_matlab::displaySearchResult_Matlab()
{
	// Call engOpen with a NULL string. This starts a MATLAB process
    // on the current host using the command "matlab".
	if (!(engine = engOpen("")))
	{
		fprintf(stderr, "\nCan't start MATLAB engine\n");
		exit(1);
	}

	// Render various Matlab plots
	displayDistributionSearchResult_Matlab();

	// Use fgetc() to make sure that we pause long enough to be
	// able to see the plot
	printf("Hit return to continue\n\n");
	fgetc(stdin);

	// Free memory
	mxDestroyArray( mxQueryDist );
	mxDestroyArray( mxNBin );


	// Close Matlab engine
	engClose( engine );

}// End function

void
FEL_matlab::displayDistributionSearchResult_Matlab( int nBin,
													int nReturnedRange,
													double* queryDist,
													double* resultRegionLimitList
													)
{
	// Create Matlab variables for our data
	mxQueryDist = mxCreateDoubleMatrix( 1, nBin, mxREAL );
	memcpy( (void *)mxGetPr( mxQueryDist ), (void *)queryDist, sizeof(double)*nBin );

	mxResultDist = mxCreateDoubleMatrix( nReturnedRange, nBin, mxREAL );
	memcpy( (void *)mxGetPr( mxResultDist ), (void *)resultDist, sizeof(double)*(nBin*nReturnedRange) );

	mxResultRegion = mxCreateDoubleMatrix( 7, nReturnedRange, mxREAL );
	memcpy( (void *)mxGetPr(mxResultRegion), (void *)resultRegionLimitList, sizeof(double)*(nReturnedRange*7) );

	double tmp = nBin;
	mxNBin = mxCreateDoubleMatrix( 1, 1, mxREAL );
	memcpy( (void *)mxGetPr( mxNBin ), (void *)(&tmp), sizeof(double) );

	double tmp2 = nReturnedRange;
	mxNReturnedRanges = mxCreateDoubleMatrix( 1, 1, mxREAL );
	memcpy( (void *)mxGetPr( mxNReturnedRanges ), (void *)(&tmp2), sizeof(double) );

	// Place the variables into the MATLAB workspace
	engPutVariable( engine, "mxQueryDist", mxQueryDist );
	engPutVariable( engine, "mxResultDist", mxResultDist );
	engPutVariable( engine, "mxNBin", mxNBin );
	engPutVariable( engine, "mxResultRegion", mxResultRegion );
	engPutVariable( engine, "mxNReturnedRanges", mxNReturnedRanges );

	// Plot the error for each domain
	engEvalString( engine, "figure" );
	engEvalString( engine, "plot( mxQueryDist, 'r' );" );
	engEvalString( engine, "title('Query-by-distribution');" );
	engEvalString( engine, "xlabel('Bin IDs');" );
	engEvalString( engine, "ylabel('Normalized Frequency');" );
	engEvalString( engine, "axis( [0 mxNBin-1 0 1] );" );

	engEvalString( engine, "figure" );
	engEvalString( engine, "mxResultDist2 = mxResultDist';" );
	engEvalString( engine, "plot( mxResultDist2(1:mxNReturnedRanges,:), 'b' );" );
	engEvalString( engine, "title('Query-by-distribution Results');" );
	engEvalString( engine, "xlabel('Bin IDs');" );
	engEvalString( engine, "ylabel('Normalized Frequency');" );
	engEvalString( engine, "axis( [0 mxNBin-1 0 1] );" );


	engEvalString( engine, "figure" );
	//engEvalString( engine, "[ nX, nY, nZ, field ] = read3DVectorField( '/media/Data/Data/FlowData/2D/2_up_2x.vec', 1, 0 );" );
	//engEvalString( engine, "[ nX, nY, nZ, field ] = read3DVectorField( '/media/Data/Data/FlowData/2D/5_up_2x.vec', 1, 0 );" );
	engEvalString( engine, "[ nX, nY, nZ, field ] = read3DVectorField( '/media/Data/Data/FlowData/2D/4_up_2x.vec', 1, 0 );" );
	engEvalString( engine, "visualize2DStreamLinesWithin2( mxResultRegion', -1, field, [nX nY nZ], [1 nX 1 nY], 20, 1, (0:1:mxNReturnedRanges) );" );

}// End function
*/



