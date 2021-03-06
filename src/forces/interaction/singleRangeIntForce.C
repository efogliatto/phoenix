#include <singleRangeIntForce.H>

using namespace std;


/** Constructor */

singleRangeIntForce::singleRangeIntForce( const string& dictName,
					  const string& eqName,
					  const latticeMesh& mesh,
					  timeOptions& Time )

    : interactionForce(dictName, eqName, mesh, Time) {



    // Update pseudo pot weights

    string model = mesh.lmodel()->name();

    if( model == "D2Q9" ) {

    	_weights.push_back( 0 );
    	_weights.push_back( 1.0/3 );
    	_weights.push_back( 1.0/3 );
    	_weights.push_back( 1.0/3 );
    	_weights.push_back( 1.0/3 );
    	_weights.push_back( 1.0/12 );
    	_weights.push_back( 1.0/12 );
    	_weights.push_back( 1.0/12 );
    	_weights.push_back( 1.0/12 );

    }

    else {

	if( model == "D3Q15" ) {

	    _weights.push_back( 0 );

	    for( uint i = 1 ; i <= 6 ; i++ )
		_weights.push_back( 1.0 / 3.0 );

	    for( uint i = 7 ; i < 15 ; i++ )
		_weights.push_back( 1.0 / 24.0 );
	    


	}

	else {

	    cout << " [ERROR]  Potential weights not yet implemented for " << model << endl;

	}

    }

}




/** Destructor */

singleRangeIntForce::~singleRangeIntForce() {}





/** Update force field. Use of virtual nodes if requested */

void singleRangeIntForce::update( scalarField& rho, scalarField& T ) {

    
    
    // Reference to neighbour array

    const vector< vector<int> >& nb = _mesh.nbArray();


    // Lattice model properties

    vector< vector<int> > vel = _mesh.lmodel()->lvel();

    vector<uint> reverse = _mesh.lmodel()->reverse();

    scalar cs2 = _mesh.lmodel()->cs2();

    scalar q = _mesh.lmodel()->q();


    

    // Move over points

    for( uint i = 0 ; i < _mesh.local() ; i++ ) {


	// Compute only if node is not on boundary

	if( _mesh.isOnBoundary(i) ) {	
    
	    for( uint j = 0 ; j < 3 ; j++ )
		_force[i][j] = 0;
		
	}
	    


	// Bulk fluid node
	
	else {

	    
	    vector<scalar> F = {0, 0, 0};	    

	    for( uint k = 1 ; k < q ; k++ ) {
       
		int neighId = nb[i][ reverse[k] ];

		scalar _rho = rho[neighId];

		scalar _T = T[neighId];
	    
		scalar alpha = _weights[k] * potential( _rho, _T, cs2 );

	    
		for( uint j = 0 ; j < 3 ; j++ ) {

		    F[j] +=  alpha * (scalar)vel[k][j] ;

		}
    

	    }

		

	    // Extra constant
		
	    scalar beta = -_G * potential( rho[i], T[i], cs2 );
    
	    for( uint j = 0 ; j < 3 ; j++) {
	
		_force[i][j] =  F[j] * beta;
	
	    }	

	    

	}
	

    }




    // Sync across processors

    _force.sync();

}
