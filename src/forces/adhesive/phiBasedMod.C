#include <phiBasedMod.H>

using namespace std;


/** Constructor */

phiBasedMod::phiBasedMod( const string& dictName,
			  const string& eqName,
			  const latticeMesh& mesh,
			  const interactionForce* Fi,
			  timeOptions& Time )

    : adhesiveForce(dictName, eqName, mesh, Fi, Time) {}



/** Destructor */

phiBasedMod::~phiBasedMod() {}



/** Update force field */

void phiBasedMod::update( scalarField& rho, scalarField& T ) {

    
    // Reference to neighbour array

    const vector< vector<int> >& nb = _mesh.nbArray();
   

    // Lattice model properties

    const vector< vector<int> >& vel = _mesh.lmodel()->lvel();

    const vector<uint>& reverse = _mesh.lmodel()->reverse();

    const scalar cs2 = _mesh.lmodel()->cs2();

    const uint q = _mesh.lmodel()->q();

    const vector<scalar>& weights = _Fi->weights();



    // Move over points

    for( uint i = 0 ; i < _mesh.local() ; i++ ) {


	// Update only for available nodes

	if( _Gads.find(i) != _Gads.end() ) {

	    
	    scalar Fads[3] = {0, 0, 0};
 
	    scalar phi = _Fi->potential( rho.at(i), T.at(i), cs2 );
    

	    for( uint k = 0 ; k < q ; k++ ) {
       
	    	int neighId = nb[i][ reverse[k] ];

	    	if(neighId == -1) {
	    
	    	    for( uint j = 0 ; j < 3 ; j++ ) {

	    		Fads[j] +=  weights[k] * (scalar)vel[k][j] ;

	    	    }

	    	}
    

	    }

    

	    // Extra constant		
    
	    for( uint j = 0 ; j < 3 ; j++)
	    	_force[i][j] =  -_Gads[i] * phi * phi * Fads[j] ;   
	    

	}


	else {

	    for( uint j = 0 ; j < 3 ; j++)
	    	_force[i][j] = 0;

	}
	

    }

   


}
