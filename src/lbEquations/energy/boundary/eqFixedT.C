#include <eqFixedT.H>

#include <random>

using namespace std;



/** Constructor */

eqFixedT::eqFixedT( const std::string& eqName,
		    const std::string& bdName,
		    const latticeMesh& mesh,	  
		    const scalarField& rho,
		    const scalarField& T,
		    const vectorField& U,
		    pdfField& pdf )
    
    : energyFixedT(eqName, bdName, mesh, rho, T, U, pdf) {}




/** Destructor */

eqFixedT::~eqFixedT() {}



/** Update pdf field */

void eqFixedT::update( const energyEquation* eeq ) {



    // Lattice constants
    
    const uint q = _mesh.lmodel()->q();

    vector<scalar> f_eq_bnd(q);

    vector<scalar> Uw = {0,0,0};




    // Move over boundary elements

    for( uint i = 0 ; i < _nodes.size() ; i++ ) {
	     	

	// Equilibrium

    	uint id = _nodes[i];

	scalar Tw = _bndVal[i];
	
	eeq->eqPS( f_eq_bnd, Tw, Uw, 0 );	


	
	// Update distribution

    	for( uint k = 0 ; k < q ; k++ ) {	    	    		       		   		    
			
	    _pdf.set(id, k, f_eq_bnd[k] );

    	}
	

    }

    

}
