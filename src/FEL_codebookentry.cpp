/*
 * FEL_codebookentry.cpp
 *
 *  Created on: Mar 8, 2012
 *      Author: abon
 */

#include "FEL_codebookentry.h"

FEL_codebookentry::FEL_codebookentry()
{
}// end constructor

FEL_codebookentry::FEL_codebookentry( int rId, double me, int nr, bool ir )
{
	rangeID = rId;
	nRot = nr;
	isRef = ir;
	matchError = me;

	rangeSegmentID = -1;
	domainSegmentID = -1;
}// end constructor

FEL_codebookentry::FEL_codebookentry( int rId, int rsId, int dId, double me, int nr, bool ir )
{
	rangeID = rId;
	rangeSegmentID = rsId;
	domainSegmentID = dId;
	nRot = nr;
	isRef = ir;
	matchError = me;

}// end constructor

FEL_codebookentry::FEL_codebookentry( const FEL_codebookentry& that )
{
	this->rangeID = that.rangeID;
	this->nRot = that.nRot;
	this->isRef = that.isRef;
	this->matchError = that.matchError;
	this->rangeSegmentID = that.rangeSegmentID;
	this->domainSegmentID = that.domainSegmentID;

}// end constructor

FEL_codebookentry&
FEL_codebookentry::operator= ( const FEL_codebookentry& that )
{
	if (this != &that ) // protect against invalid self-assignment
	{
		this->rangeID = that.rangeID;
		this->nRot = that.nRot;
		this->isRef = that.isRef;
		this->matchError = that.matchError;
		this->rangeSegmentID = that.rangeSegmentID;
		this->domainSegmentID = that.domainSegmentID;
	}
	// by convention, always return *this
	return *this;

}// end constructor

FEL_codebookentry::~FEL_codebookentry()
{
}// end destructor

void
FEL_codebookentry::setAll( int rId, double me, int nr, bool ir )
{
	rangeID = rId;
	nRot = nr;
	isRef = ir;
	matchError = me;

}// end function

int
FEL_codebookentry::getAll( int* nr, bool* ir, double* me )
{
	(*nr) = nRot;
	(*ir) = isRef;
	(*me) = matchError;

	return rangeID;

}// end function

int
FEL_codebookentry::getRangeID()
{
	return rangeID;
}// end function

int
FEL_codebookentry::getRotation()
{
	return nRot;
}// end function

bool
FEL_codebookentry::getReflection()
{
	return isRef;
}// end function


float
FEL_codebookentry::getMatchError()
{
	return matchError;

}// end function

bool
FEL_codebookentry::compare_matcherror( FEL_codebookentry cbe1, FEL_codebookentry cbe2 )
{
	double e1 = cbe1.getMatchError();
	double e2 = cbe2.getMatchError();

	if( e1 <= e2 ) return true;
	if( e1 > e2 ) return false;

}// end function
