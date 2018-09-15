#include <ppFixedU.H>

using namespace std;



/** Constructor */

ppFixedU::ppFixedU( const string& eqName, const string& bdName, const std::vector<uint>& nodes )
    : ppBndCond(nodes) {

    
    // Load boundary value

    dictionary dict("start/boundaries");

    vector<scalar> val = dict.lookUp< vector<scalar> >( eqName + "/" + bdName + "/value" );


    // Resize values at boundary and assign

    _bndVal.resize( nodes.size() );

    for( uint i = 0 ; i < _bndVal.size() ; i++ ) {

	_bndVal[i].resize(3);

	for(uint j = 0 ; j < 3 ; j++)
	    _bndVal[i][j] = val[j];

    }

    

}


/** Destructor */

ppFixedU::ppFixedU() {}









// Auxiliary eq. function

void auxEq(const latticeMesh& mesh, const scalar rho, const scalar U[3], vector<scalar>& f) {

    
    // Lattice model properties

    vector< vector<int> > vel = mesh.lmodel()->lvel();

    vector<scalar> omega = mesh.lmodel()->omega();
    
    scalar cs2 = mesh.lmodel()->cs2();

    const uint q = mesh.lmodel()->q();

    
    for( uint k = 0 ; k < q ; k++ ) {

	scalar alpha = 0,
	    beta = 0;


	// Dot product
	
	for( uint j = 0 ; j < 3 ; j++ ) {

	    alpha += vel[k][j] * U[j];

	    beta += U[j] * U[j];

	}

	    
	f[k] = rho * omega[k] * ( 1 + alpha/cs2   +   0.5 * alpha * alpha / (cs2*cs2)  -  0.5 * beta / cs2 );

    }


}





/** Update pdf field */

void ppFixedU::update( const latticeMesh& mesh, pdfField& pdf, const scalarField& rho, const scalarField& T, const vectorField& U, const pseudoPotForce* F ) {

    
    // Lattice constants
    
    const scalar q = mesh.lmodel()->q();

    vector<scalar> f_eq_nb(q);

    vector<scalar> f_eq_bnd(q);

    const vector< vector<int> >& nb = mesh.nbArray();

    const vector<uint> reverse = mesh.lmodel()->reverse();

    const vector< vector<int> > vel = mesh.lmodel()->lvel();
    


    // Move over boundary elements

    for( uint i = 0 ; i < _nodes.size() ; i++ ) {

	uint id = _nodes[i];

	scalar Uw[3] = { _bndVal[i][0], _bndVal[i][1], _bndVal[i][2] };


	for( uint k = 0 ; k < q ; k++ ) {


	    if ( nb[id][k] == -1 ) {

		
		// Need density and velocity at neighbour (reverse) node
		    
		int nbid = nb[id][ reverse[k] ];


		if( nbid != -1 ) {


		    // Compute local velocity at t+dt (macro fields are not updated yet)

		    // Local velocity
		    scalar lv[3] = {0,0,0};

		    // Local density
		    scalar nbrho = 0;

		    scalar lrho = 0;
			
	    
		    {


			// Interaction force
			
			scalar Ft[3];

			F->total(Ft,id);	


			// Move over velocity components
			for( uint jj = 0 ; jj < 3 ; jj++ ) {

			    // Move over model velocities
			    for(uint kk = 0 ; kk < q ; kk++) {

				lv[jj] += vel[kk][jj] * pdf[nbid][kk];
		    
			    }
	    
			}


			for(uint kk = 0 ; kk < q ; kk++) {

			    nbrho += pdf[nbid][kk];

			    lrho += pdf[id][kk];

			}
		

			// Add interaction force and divide by density
			for( uint jj = 0 ; jj < 3 ; jj++ ) {

			    lv[jj] = ( lv[jj]   +   Ft[jj] * 0.5  ) / nbrho;
	
			}

		    }

			
			
		    // Equilibrium

		    auxEq( mesh, nbrho, lv, f_eq_nb );

		    auxEq( mesh, nbrho, Uw, f_eq_bnd );

			
		    

		    // Update distribution
			
		    pdf.set(id, k, f_eq_bnd[k] + (pdf[nbid][k] - f_eq_nb[k] ) );

		}
		

	    }
	    

	}

    }
    
}
