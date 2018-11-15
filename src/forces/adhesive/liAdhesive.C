#include <liAdhesive.H>

#include <dictionary.H>

using namespace std;


/** Constructor */

liAdhesive::liAdhesive( const string& dictName,
			const string& eqName,
			const latticeMesh& mesh )
    
    : adhesiveForce(dictName, eqName, mesh) {


    // Read constants for each boundary

    const map< string, vector<uint> >& boundary = _mesh.boundaries();

    dictionary dict("start/boundaries");

    for( const auto &bd : boundary )
	_Gads[bd.first] = dict.lookUpOrDefault<scalar>( eqName + "/" + bd.first + "/Gads", 0 );




    // Update pseudo pot weights

    string model = mesh.lmodel()->name();

    if( model == "D2Q9" ) {

    	_weights.push_back( 0 );
    	_weights.push_back( 1.0/9 );
    	_weights.push_back( 1.0/9 );
    	_weights.push_back( 1.0/9 );
    	_weights.push_back( 1.0/9 );
    	_weights.push_back( 1.0/36 );
    	_weights.push_back( 1.0/36 );
    	_weights.push_back( 1.0/36 );
    	_weights.push_back( 1.0/36 );

    }

    else {

	if( model == "D3Q15" ) {

	    _weights.push_back( 0 );

	    for( uint i = 1 ; i <= 6 ; i++ )
		_weights.push_back( 1.0 / 9.0 );

	    for( uint i = 7 ; i < 15 ; i++ )
		_weights.push_back( 1.0 / 72.0 );
	    


	}

	else {

	    cout << " [ERROR]  Potential weights not implemented yet for " << model << endl;

	}

    }    



}


/** Destructor */

liAdhesive::~liAdhesive() {}


/** Force at specific node */

const void liAdhesive::force( const uint& i, const scalar& rho, const scalar& T, scalar f[3] ) const {


    // Initialize force

    f[0] = 0;
    f[1] = 0;
    f[2] = 0;    

    
    if( _mesh.isOnBoundary(i) ) {
	

	const scalar g_ads = _Gads.at( _mesh.nodeToBoundary(i) );

       
	// Lattice model properties

	const vector< vector<int> >& nb = _mesh.nbArray();

	const vector< vector<int> >& vel = _mesh.lmodel()->lvel();

	const vector<uint>& reverse = _mesh.lmodel()->reverse();

	const uint q = _mesh.lmodel()->q();
       


	// Compute force

	scalar presum = liAdhesive::potential(rho,T,_mesh.lmodel()->cs2());

	presum = presum * presum * (-g_ads);
	
	for( uint k = 1 ; k < q ; k++ ) {
       
	    if( nb[i][ reverse[k] ] != -1 ) {
	    
		for( uint j = 0 ; j < 3 ; j++ ) {

		    f[j] +=   _weights[k] * (scalar)vel[k][j] ;

		}

	    }
    

	}


	for( uint j = 0 ; j < 3 ; j++ )
	    f[j] = f[j] * presum;

	

    }


}




/** Interaction potential */

const scalar liAdhesive::potential( const scalar& rho, const scalar& T, const scalar& cs2 ) const {

    scalar a = 2 * (eos->p_eos(rho,T) - rho * cs2);

    scalar b(0);

    (a >= 0)  ?	 b = sqrt(a)  :  b = sqrt(-a);
    
    return b;

}
