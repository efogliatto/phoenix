#include <InamuroFixedT.H>

#include <random>

using namespace std;



/** Constructor */

InamuroFixedT::InamuroFixedT( const std::string& eqName,
			      const std::string& bdName,
			      const latticeMesh& mesh,	  
			      const scalarField& rho,
			      const scalarField& T,
			      const vectorField& U,
			      pdfField& pdf )
    
    : energyFixedT(eqName, bdName, mesh, rho, T, U, pdf) {}




/** Destructor */

InamuroFixedT::~InamuroFixedT() {}




/** Update pdf field. Inamuro */

void InamuroFixedT::update( const energyEquation* eeq ) {



    // Lattice constants
    
    const uint q = _mesh.lmodel()->q();

    const vector< vector<int> >& nb = _mesh.nbArray();

    vector<scalar> f_eq_bnd(q);

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

	for(uint j = 0 ; j < 3 ; j++)			
	    Uw[j] = _U.at(id,j);


	

	// Equilibrium populations over boundary

	// eeq->eqPS( f_eq_bnd, Tw, Uw, 0 );
	eeq->eqPS( f_eq_bnd, Tw, Uw,  eeq->heat(id) );


	
	// Update unknowk distributions

    	for( uint k = 1 ; k < q ; k++ ) {	    	    		       		   		    
			
	    if( nb[id][k] == -1 ) {

		_pdf[id][k] = f_eq_bnd[k];

	    }

    	}




	// Correction constants

	scalar beta(0), kn(0), unk(0);

    	for( uint k = 0 ; k < q ; k++ ) {	    	    		       		   		    
			
	    if( nb[id][k] == -1 ) {

		unk += _pdf[id][k];

	    }

	    else {

		kn += _pdf[id][k];

	    }

    	}


	beta = (Tw - kn) / unk;

    	for( uint k = 0 ; k < q ; k++ ) {	    	    		       		   		    
			
	    if( nb[id][k] == -1 ) {

		_pdf[id][k] = beta * _pdf[id][k];

	    }

    	}

	
	
	
	

    }

    

}
