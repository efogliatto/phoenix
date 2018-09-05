#include <scalarVector.H>

using namespace std;



/** Default constructor */

scalarVector::scalarVector() : vector<scalar>() {}



/** Constructor with size and default value */

scalarVector::scalarVector(uint sz, uint df) {

    // Allocate memory
    
    (*this).resize( sz, df );

    update();

}




/** Constructor from vector<scalar> */

scalarVector::scalarVector( const vector<scalar>& V ) {

    
    // Allocate memory
    
    (*this).resize( V.size() );



    // Copy values

    for( uint i = 0 ; i < V.size() ; i++ ) {

	(*this)[i] = V[i];

    }


    
    // Update cached size

    _size = V.size();
    

}



/** Default destructor */

scalarVector::~scalarVector() {}



/** Update vector information */

void scalarVector::update() {

    _size = vector<scalar>::size();

}



/** Overloaded << operator */

ostream& operator<<(ostream& os, const scalarVector& V) {

    for(uint i = 0 ; i < V.sz()-1 ; i++) {

	os << V[i] << " ";

    }

    os << V.back();

    return os;
    
} 
