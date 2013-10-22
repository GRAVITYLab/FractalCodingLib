/*
 * FEL_codebookentry.h
 *
 *  Created on: Mar 8, 2012
 *      Author: abon
 */

#ifndef FEL_CODEBOOKENTRY_H_
#define FEL_CODEBOOKENTRY_H_

#include <list>

class FEL_codebookentry
{
	int rangeID;
	int rangeSegmentID;
	int domainSegmentID;
	int nRot;
	bool isRef;
	double matchError;

public:

	FEL_codebookentry();
	FEL_codebookentry( int rId, double me, int nr, bool ir );
	FEL_codebookentry( int rId, int rsId, int dId, double me, int nr, bool ir );
	FEL_codebookentry( const FEL_codebookentry& that );

	FEL_codebookentry& operator= ( const FEL_codebookentry& that );

	~FEL_codebookentry();

	void setAll( int rId, double me, int nr, bool ir );

	int getRangeID();
	int getRotation();
	bool getReflection();
	float getMatchError();
	int getAll( int* nr, bool* ir, double* me );


	static bool compare_matcherror( FEL_codebookentry cbe1, FEL_codebookentry cbe2 );

};


#endif /* FEL_CODEBOOKENTRY_H_ */
