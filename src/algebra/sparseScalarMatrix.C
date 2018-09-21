#include <sparseScalarMatrix.H>

#include <algorithm>


using namespace std;


/** Default constructor */

sparseScalarMatrix::sparseScalarMatrix() {}



/** Constructor from diagonal */

sparseScalarMatrix::sparseScalarMatrix( const vector<scalar>& V ) {

    for( uint i = 0 ; i < V.size() ; i++ ) {

	addElement(V[i],i,i);

    }

}



/** Default destructor */

sparseScalarMatrix::~sparseScalarMatrix() {}
    


/** Matrix - vector multiplication */

const void sparseScalarMatrix::matDotVec (const vector<scalar>& V, vector<scalar>& res) const {

    if( res.size() == _size ) {
       	
	std::fill( res.begin(), res.end(), 0 );

	for( uint k = 0 ; k < _idx0.size() ; k++ ) {

	    uint i = _idx0[k];

	    uint j = _idx1[k];

	    res[i] += _values[k] * V[j];

	}

    }

}


    
/** Add element */

const void sparseScalarMatrix::addElement( const scalar& val, const uint& i, const uint& j ) {


    // Check if position is already set

    vector<uint>::iterator i_it = std::find(_idx0.begin(), _idx0.end(), i);

    vector<uint>::iterator j_it = std::find(_idx1.begin(), _idx1.end(), j);

    if(  ( i_it == _idx0.end() )    &&   ( j_it == _idx1.end() )  ) {

    	_values.push_back(val);

    	_idx0.push_back(i);

    	_idx1.push_back(j);

    }


    else {

    	uint k = i_it - _idx0.begin(); 

    	_values[k] = val;

    }


    
    // Update matrix size

    uint max_i = *std::max_element( _idx0.begin(), _idx0.end() );

    uint max_j = *std::max_element( _idx1.begin(), _idx1.end() );

    max_i >= max_j  ?  _size = max_i + 1  :  _size = max_j + 1;
    
}
