#include <singleRangeMixedIntForce.H>

using namespace std;


/** Constructor */

singleRangeMixedIntForce::singleRangeMixedIntForce( const string& dictName,
						    const string& eqName,
						    const latticeMesh& mesh,
						    timeOptions& Time )

    : singleRangeIntForce(dictName, eqName, mesh, Time) {



    // Read beta constant
    
    dictionary dict( dictName );

    _beta = dict.lookUp<scalar>( eqName + "/Forces/Interaction/beta" );

    

}




/** Destructor */

singleRangeMixedIntForce::~singleRangeMixedIntForce() {}



/** Update force field */

void singleRangeMixedIntForce::update( scalarField& rho, scalarField& T ) {

    
    // Reference to neighbour array

    const vector< vector<int> >& nb = _mesh.nbArray();


    // Lattice model properties

    vector< vector<int> > vel = _mesh.lmodel()->lvel();

    vector<uint> reverse = _mesh.lmodel()->reverse();

    scalar cs2 = _mesh.lmodel()->cs2();

    scalar q = _mesh.lmodel()->q();
    

    // Move over points

    for( uint i = 0 ; i < _mesh.local() ; i++ ) {


	// Check if neighbour is on boundary

	bool isOnBnd(false);

	for( uint k = 0 ; k < q ; k++ ) {

	    if( nb[i][k] == -1 ) {

		isOnBnd = true;

		k = q;

	    }

	}



	// Compute only if node is not on bounday

	if( isOnBnd ) {

	    for( uint j = 0 ; j < 3 ; j++ )
		_force[i][j] = 0;

	}


	else {

	    
	    vector<scalar> F1 = {0, 0, 0};

	    vector<scalar> F2 = {0, 0, 0};	    

	    for( uint k = 1 ; k < q ; k++ ) {
       
		int neighId = nb[i][ reverse[k] ];

		scalar _rho = rho[neighId];

		scalar _T = T[neighId];
	    
		scalar alpha_1 = _weights[k] * potential( _rho, _T, cs2 );

		scalar alpha_2 = alpha_1 * alpha_1 / _weights[k];

	    
		for( uint j = 0 ; j < 3 ; j++ ) {

		    F1[j] +=  alpha_1 * (scalar)vel[k][j] ;

		    F2[j] +=  alpha_2 * (scalar)vel[k][j] ;

		}
    

	    }

		

	    // Extra constant
		
	    scalar kappa = -_G * potential( rho[i], T[i], cs2 ) * _beta;
    
	    for( uint j = 0 ; j < 3 ; j++) {
	
		_force[i][j] =  F1[j] * kappa
		              + 0.5 * ( 1 - _beta) * F2[j] ;
	
	    }	

	    

	}
	

    }




    // Sync across processors

    _force.sync();

}
