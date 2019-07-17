#include <NEExtFixedT.H>

#include <random>

using namespace std;



/** Constructor */

NEExtFixedT::NEExtFixedT( const std::string& eqName,
			      const std::string& bdName,
			      const latticeMesh& mesh,	  
			      const scalarField& rho,
			      const scalarField& T,
			      const vectorField& U,
			      pdfField& pdf )
    
    : energyFixedT(eqName, bdName, mesh, rho, T, U, pdf) {}




/** Destructor */

NEExtFixedT::~NEExtFixedT() {}



/** Update pdf field */

void NEExtFixedT::update( const energyEquation* eeq ) {



    // Lattice constants
    
    const uint q = _mesh.lmodel()->q();

    vector<scalar> f_eq_nb(q);

    vector<scalar> f_eq_bnd(q);

    vector<scalar> Unbid = {0,0,0};

    vector<scalar> Uw = {0,0,0};


    

    // Random seeds
    
    uniform_real_distribution<scalar> unif( (100.0-_pert)/100.0, (100.0+_pert)/100.0);

    default_random_engine re;

    re.seed( _mesh.pid() );

    


    // Move over boundary elements

    for( uint i = 0 ; i < _nodes.size() ; i++ ) {
	

    	uint id = _nodes[i];

	scalar Tw = unif(re) * _bndVal[i];


	
	// Density and velocity at neighbour node

	const uint nbid = _nbid[i];
	


	// Velocity at boundary

	for(uint j = 0 ; j < 3 ; j++) {
			
	    Unbid[j] = _U.at(nbid,j);

	    // Uw[j] = _U.at(id,j);

	}

	

	// Equilibrium

	eeq->eqPS( f_eq_nb, _T.at(nbid), Unbid, 0 );

	eeq->eqPS( f_eq_bnd, Tw, Uw, 0 );

	// eeq->eqPS( f_eq_nb, _T.at(nbid), Unbid, eeq->heat(id) );

	// eeq->eqPS( f_eq_bnd, Tw, Uw, eeq->heat(id) );
	


	
	// Update distribution

    	for( uint k = 0 ; k < q ; k++ ) {	    	    		       		   		    
			
	    _pdf[id][k] = f_eq_bnd[k] + _pdf[nbid][k] - f_eq_nb[k];

    	}
	

    }

    

}
