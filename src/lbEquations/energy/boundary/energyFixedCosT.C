#include <energyFixedCosT.H>

#include <math.h>

using namespace std;


/** Constructor */

energyFixedCosT::energyFixedCosT( const std::string& eqName,
				    const std::string& bdName,
				    const latticeMesh& mesh,
				    const scalarField& rho,
				    const scalarField& T,
				    const vectorField& U,
				    pdfField& pdf )

    : energyFixedT( eqName, bdName, mesh, rho, T, U, pdf ) {
    

    
    // Oscilation parameters

    dictionary dict("start/boundaries");
    
    scalar A = dict.lookUp<scalar>( eqName + "/" + bdName + "/amplitude" );

    vector<scalar> w = dict.lookUp< vector<scalar> >( eqName + "/" + bdName + "/length" );

    for( uint i = 0 ; i < 3 ; i++ ) {
	
	if(w[i] != 0)
	    w[i] = 2*M_PI/w[i];	    

    }


    // Reassign boundary values
    
    for( uint i = 0 ; i < _bndVal.size() ; i++ ) {

	const vector<int>& lp = _mesh.latticePoint( _nodes[i] );
	
	scalar b = w[0]*lp[0] + w[1]*lp[1] + w[2]*lp[2];
	
	_bndVal[i] = _bndVal[i]*( 1 + A*cos(b) );

    }

    

    


    
}




/** Destructor */

energyFixedCosT::~energyFixedCosT() {}



/** Update pdf field */

void energyFixedCosT::update( const energyEquation* eeq ) {
    
    energyFixedT::update( eeq );

}
