#include <scalarMatrix.H>

using namespace std;



/** Default constructor */

scalarMatrix::scalarMatrix() : vector< vector<scalar> >() {}



/** Constructor from vector< vector<scalar> > */

scalarMatrix::scalarMatrix( const vector< vector<scalar> >& M ) {

    
    // Allocate memory
    
    (*this).resize( M.size() );

    for( uint i = 0 ; i < M.size() ; i++ ) {

	(*this)[i].resize( M[i].size() );

    }



    // Copy values

    for( uint i = 0 ; i < M.size() ; i++ ) {

	for( uint j = 0 ; j < M[i].size() ; j++ ) {

	    (*this)[i][j] = M[i][j];

	}

    }


    
    // Update cached size

    _size = M.size();
    

}



/** Default destructor */

scalarMatrix::~scalarMatrix() {}






/** Update matrix information */

void scalarMatrix::update() {

    _size = vector< vector<scalar> >::size();

}



/** Overloaded << operator */

ostream& operator<<(ostream& os, const scalarMatrix& M) {

    for(uint i = 0 ; i < M.sz() ; i++) {

	for(uint j = 0 ; j < M.sz() ; j++) {

	    os << M[i][j] << " ";

	}

	os << endl;

    }    

    return os;
    
} 



/** Matrix - vector multiplication */

const void scalarMatrix::matDotVec (const scalarVector& V, scalarVector& res) const {

    
    // Res MUST be preallocated for efficiency

    if( _size == V.sz() ) {

	for( uint i = 0 ; i < _size ; i++ ) {

	    res[i] = 0;

	    for( uint j = 0 ; j < _size ; j++ ) {

		res[i] += (*this)[i][j] * V[j];

	    }

	}	

    }

}
