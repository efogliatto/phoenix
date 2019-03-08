#include <ppNEExt.H>

using namespace std;



/** Constructor */

ppNEExt::ppNEExt( const std::string& eqName,
		const std::string& bdName,
		const latticeMesh& mesh,	  
		const scalarField& rho,
		const scalarField& T,
		const vectorField& U,
		pdfField& pdf )
    
    : ppWetNodeBnd(eqName, bdName, "NEExt", mesh, rho, T, U, pdf) {

    
}


/** Destructor */

ppNEExt::~ppNEExt() {}




/** Update pdf field */

void ppNEExt::update( const pseudoPotEquation* ppeq ) {

    
    // Lattice constants
    
    const uint q = _mesh.lmodel()->q();

    vector<scalar> f_eq_nb(q);

    vector<scalar> f_eq_bnd(q);

    


    // Move over boundary elements

    for( uint i = 0 ; i < _nodes.size() ; i++ ) {
	

    	uint id = _nodes[i];

    	vector<scalar> Uw = { _bndVal[i][0], _bndVal[i][1], _bndVal[i][2] };

	vector<scalar> nvel = { 0,0,0 };


	// Density and velocity at neighbour node

	const uint nbid = _nbid[i];
	
	scalar nbrho = ppeq->localDensity(nbid);	

	ppeq->localVelocity(nvel, nbid, true);



	// Density at wall

	scalar rhow = ppeq->localDensityWithUnknowns( id, _normal[i] );
	


	// Equilibrium

	ppeq->eqPS( f_eq_nb, nbrho, nvel );

	ppeq->eqPS( f_eq_bnd, rhow, Uw );



	
	// Update distribution

    	for( uint k = 0 ; k < q ; k++ ) {	    	    		       		   		    
			
	    _pdf.set(id, k, f_eq_bnd[k] + (_pdf[nbid][k] - f_eq_nb[k] ) );

    	}
	

    }
    
}
