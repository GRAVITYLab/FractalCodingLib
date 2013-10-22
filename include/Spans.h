/*
 * Spans.h
 *
 *  Created on: Mar 20, 2013
 *      Author: abon
 */

#ifndef SPANS_H_
#define SPANS_H_

#include <cstdio>
#include <iostream>
#include <cmath>

struct MapValue
{
	int size;
	int encodingIndex;

	MapValue( int s, int e )
	{
		size = s;
		encodingIndex = e;
	}
};

struct Span1D
{
	int low;
	int high;
	float center;
	int size;

	Span1D( int x, int X )
	{
		low = x; high = X;

		center = ( low + high ) / 2.0;

		size = ( high - low + 1 );
	}

	int getSize()
	{
		return (high-low+1);
	}

	static int
	hash_int( int a)
	{
	    a -= (a<<6);
	    a ^= (a>>17);
	    a -= (a<<9);
	    a ^= (a<<4);
	    a -= (a<<3);
	    a ^= (a<<10);
	    a ^= (a>>15);
	    return a;
	}

	bool operator< ( const Span1D& n ) const
	{
		if ( this->low == n.low &&
		 this->high == n.high )
			return false;

		if( this->center < n.center )
			return true;
		else if( this->center == n.center &&
				  this->size < n.size )
			return true;

	   return false;
	}
};

struct Span2D
{
	int low[2];
	int high[2];
	float center[2];
	int size;

	Span2D( int x, int X, int y, int Y )
	{
		low[0] = x; high[0] = X;
		low[1] = y; high[1] = Y;

		center[0] = ( low[0] + high[0] ) / 2.0;
		center[1] = ( low[1] + high[1] ) / 2.0;

		size = ( high[0] - low[0] + 1 )*
				( high[1] - low[1] + 1 );
	}

	int getSize()
	{
		return (high[0]-low[0]+1)*
				(high[1]-low[1]+1);
	}

	bool operator< ( const Span2D& n ) const
	{
		if ( this->low[0] == n.low[0] &&
			 this->low[1] == n.low[1] &&
			 this->high[0] == n.high[0] &&
			 this->high[1] == n.high[1] )
		   		return false;

			if( this->center[0] < n.center[0] )
				return true;
			else if( this->center[0] == n.center[0] &&
					  this->center[1] < n.center[1] )
				return true;
			else if( this->center[0] == n.center[0] &&
					  this->center[1] == n.center[1] &&
					  this->size < n.size )
				return true;

	    	return false;

	   return false;
	}

};

struct Span3D
{
	int low[3];
	int high[3];
	float center[3];
	int size;

	Span3D( int x, int X, int y, int Y, int z, int Z )
	{
		low[0] = x; high[0] = X;
		low[1] = y; high[1] = Y;
		low[2] = z; high[2] = Z;

		center[0] = ( low[0] + high[0] ) / 2.0;
		center[1] = ( low[1] + high[1] ) / 2.0;
		center[2] = ( low[2] + high[2] ) / 2.0;

		size = ( high[0] - low[0] + 1 )*
				( high[1] - low[1] + 1 )*
				( high[2] - low[2] + 1 );

	}

	int
	getSize()
	{
		return (high[0]-low[0]+1)*
				(high[1]-low[1]+1)*
				(high[2]-low[2]+1);
	}

	void
	print()
	{
		fprintf( stderr, "Span: [%d %d %d] - [%d %d %d]\n",
						low[0], low[1], low[2],
						high[0], high[1], high[2] );
	}


	bool operator< ( const Span3D& n ) const
	{
		// Special case: to make sure equality
		// cases hold in both directions
		if ( this->low[0] == n.low[0] &&
	   		 this->low[1] == n.low[1] &&
	   		 this->low[2] == n.low[2] &&
	   		 this->high[0] == n.high[0] &&
	   		 this->high[1] == n.high[1] &&
	   		 this->high[2] == n.high[2] )
	   		return false;

		//std::cout << center[0] << ":" << center[1] << ":" << center[2] << std::endl;
		//std::cout << n.center[0] << ":" << n.center[1] << ":" << n.center[2] << std::endl;


		if( this->center[0] < n.center[0] )
			return true;
		else if( this->center[0] == n.center[0] &&
				  this->center[1] < n.center[1] )
			return true;
		else if( this->center[0] == n.center[0] &&
				  this->center[1] == n.center[1] &&
				  this->center[2] < n.center[2] )
			return true;
		else if( this->center[0] == n.center[0] &&
				  this->center[1] == n.center[1] &&
				  this->center[2] == n.center[2] &&
				  this->size < n.size )
			return true;

  		return false;
	}
};

struct Span3DComp
{
  bool
  operator() ( const Span3D& lhs, const Span3D& rhs) const
  {
	  return ( lhs.low[0] < rhs.low[0]) ||
			  ( ( lhs.low[0] == rhs.low[0]) && ( lhs.low[1] < rhs.low[1]) ) ||
			  ( ( lhs.low[0] == rhs.low[0]) && ( lhs.low[1] == rhs.low[1]) && ( lhs.low[2] < rhs.low[2]) );
  }
};

struct hashSpan1D {
    size_t operator() ( const Span1D& k) const
    {
    	int a = (int)floor( ( 2*k.low + 3*k.high ) / 2.0f  );

		return (51 + Span1D::hash_int( a ) );
    }
};

struct equalSpan1D {
    bool operator() ( const Span1D& id1, const Span1D& id2 ) const
    {
		return ( id1.low == id2.low );
    }
};

struct hashSpan2D {
    size_t operator() ( const Span2D& k) const
    {
		int a = (int)floor( ( 2*k.low[0] + 3*k.high[0] ) / 2.0 );
		int b = (int)floor( ( 2*k.low[1] + 3*k.high[1] ) / 2.0 );

		return ( (51 + Span1D::hash_int( a ) * 51 + Span1D::hash_int( b ) ) );
    }
};

struct equalSpan2D {
    bool operator() ( const Span2D& id1, const Span2D& id2 ) const
    {
		return ( id1.low[0] == id2.low[0] ) &&
				( id1.low[1] == id2.low[1] ) &&
				( id1.high[0] == id2.high[0] ) &&
				( id1.high[1] == id2.high[1] );

    }
};

struct hashSpan3D {
    size_t operator() ( const Span3D& k) const
    {
		int a = (int)floor( ( 2*k.low[0] + 3*k.high[0] ) / 2.0 );
		int b = (int)floor( ( 2*k.low[1] + 3*k.high[1] ) / 2.0 );
		int c = (int)floor( ( 2*k.low[2] + 3*k.high[2] ) / 2.0 );

		return ( (51 + Span1D::hash_int( a ) * 51 + Span1D::hash_int( b ) ) + Span1D::hash_int( c )  );
    }
};

struct equalSpan3D {
    bool operator() ( const Span3D& id1, const Span3D& id2 ) const
    {
		return ( id1.low[0] == id2.low[0] ) &&
				( id1.low[1] == id2.low[1] ) &&
				( id1.low[2] == id2.low[2] ) &&
				( id1.high[0] == id2.high[0] ) &&
				( id1.high[1] == id2.high[1] ) &&
				( id1.high[2] == id2.high[2] );

    }
};

struct lessThanSpan3D {

	bool operator() ( const Span3D& left, const Span3D& right ) const
    {
		#ifdef DEBUG_MODE
		fprintf( stderr, "left:[%d %d %d][%d %d %d]\n",
				left.low[0], left.low[1], left.low[2],
				left.high[0], left.high[1], left.high[2] );
		fprintf( stderr, "right:[%d %d %d][%d %d %d]\n",
				right.low[0], right.low[1], right.low[2],
				right.high[0], right.high[1], right.high[2] );
		#endif

		/*
		if ( left.low[0] == right.low[0] &&
	   		 left.low[1] == right.low[1] &&
	   		 left.low[2] == right.low[2] &&
	   		 left.high[0] == right.high[0] &&
	   		 left.high[1] == right.high[1] &&
	   		 left.high[2] == right.high[2] )
	   		return false;
	   		*/

		if( left.center[0] < right.center[0] )
			return true;
		else if( left.center[0] == right.center[0] &&
				  left.center[1] < right.center[1] )
			return true;
		else if( left.center[0] == right.center[0] &&
				  left.center[1] == right.center[1] &&
				  left.center[2] < right.center[2] )
			return true;
		else if( left.center[0] == right.center[0] &&
				  left.center[1] == right.center[1] &&
				  left.center[2] == right.center[2] &&
				  left.size < right.size )
			return true;

    	return false;

   		/*
		if ( left.high[0] < right.low[0] )
		   return true;
   		if ( left.high[0] == right.low[0] &&
   		     left.low[0] < right.high[0] )
		   return true;
    	if ( left.low[0] <= right.low[0] &&
    		 right.high[0] < left.high[0] )
		   return true;
    	if ( left.low[0] < right.low[0] &&
    		 right.high[0] <= left.high[0] )
		   return true;
    	if ( left.low[0] < right.low[0] &&
    		 left.high[0] <= right.high[0] )
		   return true;
    	if ( left.low[0] <= right.low[0] &&
    		 left.high[0] < right.high[0] )
		   return true;

   		if ( left.high[1] < right.low[1] )
		   return true;
   		if ( left.high[1] == right.low[1] &&
   		     left.low[1] < right.high[1] )
		   return true;
    	if ( left.low[1] <= right.low[1] &&
    		 right.high[1] < left.high[1] )
		   return true;
    	if ( left.low[1] < right.low[1] &&
    		 right.high[1] <= left.high[1] )
		   return true;
    	if ( left.low[1] < right.low[1] &&
    		 left.high[1] <= right.high[1] )
		   return true;
    	if ( left.low[1] <= right.low[1] &&
    		 left.high[1] < right.high[1] )
		   return true;

   		if ( left.high[2] < right.low[2] )
		   return true;
   		if ( left.high[2] == right.low[2] &&
   		     left.low[2] < right.high[2] )
		   return true;
    	if ( left.low[2] <= right.low[2] &&
    		 right.high[2] < left.high[2] )
		   return true;
    	if ( left.low[2] < right.low[2] &&
    		 right.high[2] <= left.high[2] )
		   return true;
    	if ( left.low[2] < right.low[2] &&
    		 left.high[2] <= right.high[2] )
		   return true;
    	if ( left.low[2] <= right.low[2] &&
    		 left.high[2] < right.high[2] )
		   return true;
		   */

    }// end operator
};


#endif /* SPANS_H_ */
