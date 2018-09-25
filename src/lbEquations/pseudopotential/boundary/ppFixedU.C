#include <ppFixedU.H>

using namespace std;



/** Constructor */

ppFixedU::ppFixedU( const std::string& eqName,
		    const std::string& bdName,
		    const latticeMesh& mesh,	  
		    const scalarField& rho,
		    const scalarField& T,
		    const vectorField& U,
		    pdfField& pdf )
    
    : ppBndCond(mesh, rho, T, U, pdf, bdName, "fixedU") {

    
    // Load boundary value

    dictionary dict("start/boundaries");

    vector<scalar> val = dict.lookUp< vector<scalar> >( eqName + "/" + bdName + "/value" );


    // Resize values at boundary and assign

    _bndVal.resize( _nodes.size() );

    for( uint i = 0 ; i < _bndVal.size() ; i++ ) {

    	_bndVal[i].resize(3);

    	for(uint j = 0 ; j < 3 ; j++)
    	    _bndVal[i][j] = val[j];

    }

    

}


/** Destructor */

ppFixedU::~ppFixedU() {}





/** Update pdf field */

void ppFixedU::update( const pseudoPotEquation* ppeq ) {

    
    // Lattice constants
    
    const uint q = _mesh.lmodel()->q();

    vector<scalar> f_eq_nb(q);

    vector<scalar> f_eq_bnd(q);

    const vector< vector<int> >& nb = _mesh.nbArray();

    const vector<uint> reverse = _mesh.lmodel()->reverse();

    const vector< vector<int> > vel = _mesh.lmodel()->lvel();
    


    // Move over boundary elements

    for( uint i = 0 ; i < _nodes.size() ; i++ ) {
	

    	uint id = _nodes[i];

    	vector<scalar> Uw = { _bndVal[i][0], _bndVal[i][1], _bndVal[i][2] };

	vector<scalar> nvel = { 0,0,0 };


    	for( uint k = 0 ; k < q ; k++ ) {


    	    if ( nb[id][k] == -1 ) {

		
    		// Need density and velocity at neighbour (reverse) node
		    
    		int nbid = nb[id][ reverse[k] ];


    		if( nbid != -1 ) {


    		    // Density and velocity at neighbour node
		    
    		    scalar nbrho = ppeq->localDensity(nbid);

		    ppeq->localVelocity(nvel, nbid, true);

	    
	    			
    		    // Equilibrium

    		    ppeq->eqPS( f_eq_nb, nbrho, nvel );

    		    ppeq->eqPS( f_eq_bnd, nbrho, Uw );
		    
		    

    		    // Update distribution
			
    		    _pdf.set(id, k, f_eq_bnd[k] + (_pdf[nbid][k] - f_eq_nb[k] ) );

    		}
		

    	    }
	    

    	}

    }
    
}
