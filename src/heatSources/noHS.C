#include <noHS.H>

using namespace std;


/** Constructor */

noHS::noHS( const string& dictName,
	    const string& eqName,
	    const latticeMesh& mesh,
	    timeOptions& Time )

    : heatSource(dictName, eqName, mesh, Time) {       


    // Move over points

    for( uint id = 0 ; id < _mesh.npoints() ; id++ ) {

    	_source[id] = 0;       
    
    }
    
}




/** Update source field */

void noHS::update( const scalarField& rho, const scalarField& T, const vectorField& U ) {}
