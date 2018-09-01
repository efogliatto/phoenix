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



/** Matrix - vector multiplication */

const scalarVector scalarMatrix::operator* (const scalarVector& V) const{

    scalarVector result( V.sz(), 0 );

    // result.resize(V.sz());
    
    // for( uint i = 0 ; i < _size ; i++) {
    // 	result[i] = 0;
    // }

    
    if( _size == V.sz() ) {

	for( uint i = 0 ; i < _size ; i++) {

	    for( uint j = 0 ; j < _size ; j++) {

		result[i] += (*this)[i][j] * V[j];

	    }

	}

    }


    result.update();

    return result;

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
