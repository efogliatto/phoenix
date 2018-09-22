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



    dictionary prop("properties/macroProperties");

    _a1 = prop.lookUp<scalar>( eqName + "/HeatSource/Constants/alpha_1");

    _a2 = prop.lookUp<scalar>( eqName + "/HeatSource/Constants/alpha_2");

}



/** Destructor */

energyFixedT::~energyFixedT() {}




void eq( vector<scalar>& m, const scalar& T_, const vector<scalar>& U_, const latticeMesh& mesh, const scalar& a1, const scalar& a2 ) {

    
    const scalar _U[3] = { U_[0], U_[1], U_[2] };

    const uint q = mesh.lmodel()->q();

    vector<scalar> n_eq(q);

    const scalarMatrix& invM = mesh.lmodel()->MRTInvMatrix();
    

    
    switch(q) {

    case 9:

	n_eq[0] = T_;
	
    	n_eq[1] = a1 * T_;
	
    	n_eq[2] = a2 * T_;
	
    	n_eq[3] = T_ * _U[0];
	
    	n_eq[4] = T_ * (-_U[0]);
	
    	n_eq[5] = T_ * _U[1];
	
    	n_eq[6] = T_ * (-_U[1]);
	
    	n_eq[7] = 0;
	
    	n_eq[8] = 0;


    	break;


    case 15:

	n_eq[0]  = T_;
	
	n_eq[1]  = a1 * T_;

	n_eq[2]  = a2 * T_;

	n_eq[3]  =   T_ * _U[0];

	n_eq[4]  = -(T_ * _U[0]);

	n_eq[5]  =   T_ * _U[1];

	n_eq[6]  = -(T_ * _U[1]);

	n_eq[7]  =   T_ * _U[2];

	n_eq[8]  = -(T_ * _U[2]);

	n_eq[9]  = 0;

	n_eq[10] = 0;

	n_eq[11] = 0;

	n_eq[12] = 0;

	n_eq[13] = 0;

	n_eq[14] = 0;
	
    	break;

    }



    invM.matDotVec(n_eq,m);


}



/** Update pdf field */

void energyFixedT::update( const latticeMesh& mesh ) {


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
		    

		    eq( f_eq_nb, _T.at(nbid), Unbid, mesh, _a1, _a2 );

		    eq( f_eq_bnd, Tw, Uw, mesh, _a1, _a2 );
	    

		    
		    // Update distribution
			
		    _pdf.set(id, k, f_eq_bnd[k] + (_pdf[nbid][k] - f_eq_nb[k] ) );		    

		}
		

	    }
	    

	}


    }

}
