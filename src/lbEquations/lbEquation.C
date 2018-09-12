#include <lbEquation.H>

using namespace std;


/** Default constructor */

lbEquation::lbEquation( const string& name,
			const latticeMesh& mesh_,
			timeOptions& Time_,
			pdfField& pdf_) : ename (name),
					  mesh(mesh_),
					  Time(Time_),
					  _q( mesh.lmodel()->q() ),					  
					  pdf(pdf_) {


    // Resize swap array

    swap.resize( mesh.local() );

    for( uint i = 0 ; i < swap.size() ; i++ )
	swap[i].resize(_q);
    

}




/** Default destructor */

lbEquation::~lbEquation() {}



/** Collision process */

const void lbEquation::collision() {}




/** Streamming process */

const void lbEquation::streaming() {


    // References

    const vector< vector<int> > nb = mesh.nbArray();
    
    
    
    // Copy all values to swap

    for( uint i = 0 ; i < mesh.local() ; i++ ) {

	for( uint k = 0 ; k < _q ; k++ ) {

	    swap[i][k] = pdf[i][k];
	    
	}

    }


    // Copy only neighbours to swap   

    for( uint i = 0 ; i < mesh.local() ; i++ ) {

    	for( uint k = 0 ; k < _q ; k++ ) {

    	    int neighId = nb[i][k];

    	    if( neighId != -1 ) {

    		swap[i][k] = pdf[neighId][k];

    	    }
    
    	}

    }



    // Copy back from swap
    
    for( uint i = 0 ; i < mesh.local() ; i++ ) {

    	for( uint k = 0 ; k < _q ; k++ ) {

    	    pdf[i][k] = swap[i][k];
	    
    	}

    }


}



/** Set pdf to equilibrium values */

const void lbEquation::setEquilibrium() {}
