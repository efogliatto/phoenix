#include <liAdhesive.H>

#include <dictionary.H>

using namespace std;


/** Constructor */

liAdhesive::liAdhesive( const string& dictName,
			const string& eqName,
			const latticeMesh& mesh )
    
    : simpleAdhesive(dictName, eqName, mesh) {}


/** Destructor */

liAdhesive::~liAdhesive() {}


/** Force at specific node */

const void liAdhesive::force( const uint& i, const scalarField& rho, const scalarField& T, scalar f[3] ) const {

    f[0] = 0;
    f[1] = 0;
    f[2] = 0;
    
    
    if( closestNodes.find(i) != closestNodes.end() ) {

	
    	// Reference to neighbour array

    	const vector< vector<int> >& nb = _mesh.nbArray();


    	// Lattice model properties

    	const vector< vector<int> >& vel = _mesh.lmodel()->lvel();

    	const vector<uint>& reverse = _mesh.lmodel()->reverse();

    	const scalar cs2 = _mesh.lmodel()->cs2();

    	const uint q = _mesh.lmodel()->q();


    	vector<scalar> F = {0, 0, 0};	    

    	for( uint k = 1 ; k < q ; k++ ) {
       
    	    int neighId = nb[i][ reverse[k] ];

	    if( _mesh.isOnBoundary(neighId) ) {
	    
		for( uint j = 0 ; j < 3 ; j++ ) {

		    F[j] +=  _weights[k] * (scalar)vel[k][j] ;

		}

	    }
    

    	}

		

	// Extra constant
		
	scalar beta = potential( rho.at(i), T.at(i), cs2 );

	beta = beta * closestNodes.at(i);	    
    
	for( uint j = 0 ; j < 3 ; j++)
	    f[j] =  F[j] * beta;
	
	

    }

}
