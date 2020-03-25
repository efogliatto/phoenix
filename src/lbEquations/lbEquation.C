#include <lbEquation.H>

using namespace std;


/** Default constructor */

lbEquation::lbEquation( const string& name,
			const latticeMesh& mesh_,
			timeOptions& Time_,
			pdfField& pdf_)
    : ename (name),
      mesh(mesh_),
      Time(Time_),
      _pdf(pdf_) {



    // Read coefficients from dictionary (DEPRECATED)

    #pragma message( "Single Tau use is DEPRECATED" )

    dictionary dict("properties/macroProperties");

    _Tau = dict.lookUp< vector<scalar> >( name + "/LBModel/Tau" );


    

    // Resize swap array

    const uint q = mesh.lmodel()->q();

    const uint np = mesh.local();

    _swap.resize( np );

    for( uint i = 0 ; i < np ; i++ )
    	_swap[i].resize(q);



    // Relaxation model

    relaxModelCreator rmodel;

    _relax = rmodel.create( name + "/LBModel" );
    
    
}




/** Default destructor */

lbEquation::~lbEquation() {}



/** Collision process */

const void lbEquation::collision() {}




/** Streamming process */

const void lbEquation::streaming() {


    // References

    const vector< vector<int> >& nb = mesh.nbArray();

    const uint q = mesh.lmodel()->q();
    
    

    // Copy only neighbours to swap   

    for( uint i = 0 ; i < mesh.local() ; i++ ) {

    	for( uint k = 0 ; k < q ; k++ ) {

    	    int neighId = nb[i][k];

    	    if( neighId != -1 ) {

    		_swap[i][k] = _pdf[neighId][k];

    	    }

	    else {

		_swap[i][k] = _pdf[i][k];

	    }
    
    	}

    }



    // Copy back from swap
    
    for( uint i = 0 ; i < mesh.local() ; i++ ) {

    	for( uint k = 0 ; k < q ; k++ ) {

    	    _pdf[i][k] = _swap[i][k];
	    
    	}

    }


}



/** Set pdf to equilibrium values */

const void lbEquation::setEquilibrium() {}
