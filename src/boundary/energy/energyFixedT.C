#include <energyFixedT.H>

using namespace std;



/** Constructor */

energyFixedT::energyFixedT( const std::string& eqName,
			    const std::string& bdName,
			    const std::vector<uint>& nodes,
			    const scalarField& rho,
			    const scalarField& T,
			    const vectorField& U,
			    pdfField& pdf )
    
    : energyBndCond(nodes, rho, T, U, pdf) {

    

    // Load boundary value

    dictionary dict("start/boundaries");

    scalar val = dict.lookUp<scalar>( eqName + "/" + bdName + "/value" );


    // Resize values at boundary and assign

    _bndVal.resize( nodes.size() );

    for( uint i = 0 ; i < _bndVal.size() ; i++ )
	    _bndVal[i] = val;

    

}



/** Destructor */

energyFixedT::~energyFixedT() {}



/** Update pdf field */

void energyFixedT::update( const latticeMesh& mesh, std::function<void(std::vector<scalar>&, const scalar&, const vector<scalar>&)> eq ) {

    
    // Lattice constants
    
    const scalar q = mesh.lmodel()->q();

    vector<scalar> f_eq_nb(q);

    vector<scalar> f_eq_bnd(q);

    const vector< vector<int> >& nb = mesh.nbArray();

    const vector<uint> reverse = mesh.lmodel()->reverse();

    vector<scalar> Unbid = {0,0,0};

    vector<scalar> Uw = {0,0,0};

    


    // Move over boundary elements

    for( uint i = 0 ; i < _nodes.size() ; i++ ) {

	uint id = _nodes[i];

	scalar Tw = _bndVal[i];


	for( uint k = 0 ; k < q ; k++ ) {


	    if ( nb[id][k] == -1 ) {

		
		// Need density and velocity at neighbour (reverse) node
		    
		int nbid = nb[id][ reverse[k] ];


		if( nbid != -1 ) {

			
		    // Equilibrium

		    for(uint j = 0 ; j < 3 ; j++) {
			
			Unbid[j] = _U.at(nbid,j);

			Uw[j] = _U.at(id,j);

		    }
		    

		    eq( f_eq_nb, _T.at(nbid), Unbid );

		    eq( f_eq_bnd, Tw, Uw );
	    

		    
		    // Update distribution
			
		    _pdf.set(id, k, f_eq_bnd[k] + (_pdf[nbid][k] - f_eq_nb[k] ) );		    

		}
		

	    }
	    

	}


    }

}
