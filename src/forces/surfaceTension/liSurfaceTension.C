#include <liSurfaceTension.H>

#include <dictionary.H>

using namespace std;


/** Constructor */

liSurfaceTension::liSurfaceTension( const string& dictName,
				    const string& eqName,
				    const latticeMesh& mesh )
    
    : surfaceTension(dictName, eqName, mesh), _kappa(0) {


    
    // Update model constant
    
    dictionary dict(dictName);

    _kappa = dict.lookUp<scalar>( eqName + "/SurfaceTension/kappa" );
    

}




/** Destructor */

liSurfaceTension::~liSurfaceTension() {}






/** Force at specific node */

void liSurfaceTension::ST( const uint& i, const scalarField& rho, const scalarField& T, vector<scalar>& C, interactionForce* _fi, const vector<scalar>& Tau ) {


    // Lattice constants
    
    const uint q = _mesh.lmodel()->q();

    const scalar cs2 = _mesh.lmodel()->cs2();

    const vector< vector<int> > vel = _mesh.lmodel()->lvel();

    const vector<uint>& reverse = _mesh.lmodel()->reverse();

    const vector< vector<int> >& nb = _mesh.nbArray();


    // Initialize term
    
    for(uint j = 0 ; j < q ; j++)
	C[j] = 0;



   
    
    if( q == 9 ) {


	// Do not use unexisting neighbour
	if( _mesh.isOnBoundary(i) == false ) {


	    // Q tensor

	    scalar Q[2][2] = {{0,0},{0,0}};

	    const scalar phi = _fi->potential(rho.at(i), T.at(i), cs2);

	    const scalar omega[9] = { 0, 1.0/3, 1.0/3, 1.0/3, 1.0/3, 1.0/12, 1.0/12, 1.0/12, 1.0/12 };
	


	    // Move over velocities

	    for( uint k = 1 ; k < q ; k++ ) {

	    	int neighId = nb[i][ reverse[k] ];

		scalar alpha(0);
		
	    	if (neighId != -1) {

	    	    alpha = omega[k] * ( _fi->potential( rho.at(neighId), T.at(neighId), cs2 ) - phi );

	    	}


	    	for( uint ii = 0 ; ii < 2 ; ii++ ) {

	    	    for( uint jj = 0 ; jj < 2 ; jj++ ) {

	    		Q[ii][jj] += alpha * vel[k][ii] * vel[k][jj];

	    	    }

	    	}
		    

	    }


	    for( uint ii = 0 ; ii < 2 ; ii++ ) {

	    	for( uint jj = 0 ; jj < 2 ; jj++ ) {

	    	    Q[ii][jj] = -Q[ii][jj] * 0.5 * phi * _kappa;

	    	}

	    }




	    // Complete C array

	    C[0] = 0;

	    C[1] = 1.5 * Tau[1] * (Q[0][0] + Q[1][1]);

	    C[2] = -1.5 * Tau[2] * (Q[0][0] + Q[1][1]);

	    C[3] = 0;

	    C[4] = 0;

	    C[5] = 0;

	    C[6] = 0;

	    C[7] = -Tau[8] * (Q[0][0] - Q[1][1]);

	    C[8] = -Tau[8] * Q[0][1];	    



	}


	
	// Boundary node
	    
	else {

	    for( uint k = 0 ; k < q ; k++ ) {

		C[k] = 0;

	    }

	}




	

    }

    else {

	if( q == 15 ) {




	    // Do not use unexisting neighbour
	    if( _mesh.isOnBoundary(i) == false ) {


		// Q tensor

		scalar Q[3][3] = { {0,0,0}, {0,0,0}, {0,0,0} };

		const scalar phi = _fi->potential(rho.at(i), T.at(i), cs2);

		const scalar omega[15] = { 0, 1.0/3, 1.0/3, 1.0/3, 1.0/3, 1.0/3, 1.0/3, 1.0/24, 1.0/24, 1.0/24, 1.0/24, 1.0/24, 1.0/24, 1.0/24, 1.0/24  };
	


		// Move over velocities

		for( uint k = 1 ; k < q ; k++ ) {

		    int neighId = nb[i][ reverse[k] ];

		    scalar alpha(0);
		
		    if (neighId != -1) {

			alpha = omega[k] * ( _fi->potential( rho.at(neighId), T.at(neighId), cs2 ) - phi );

		    }


		    for( uint ii = 0 ; ii < 3 ; ii++ ) {

			for( uint jj = 0 ; jj < 3 ; jj++ ) {

			    Q[ii][jj] += alpha * vel[k][ii] * vel[k][jj];

			}

		    }
		    

		}


		for( uint ii = 0 ; ii < 3 ; ii++ ) {

		    for( uint jj = 0 ; jj < 3 ; jj++ ) {

			Q[ii][jj] = -Q[ii][jj] * 0.5 * phi * _kappa;

		    }

		}




		// Complete C array

		C[0] = 0;

		C[1] = (4.0/5.0) * Tau[1] * (Q[0][0] + Q[1][1] + Q[2][2]);

		C[2] = 0;

		C[3] = 0;

		C[4] = 0;

		C[5] = 0;

		C[6] = 0;

		C[7] = 0;

		C[8] = 0;

		C[9] = -Tau[10]*( 2.0*Q[0][0] - Q[1][1] - Q[2][2] );

		C[10] = -Tau[10]*( Q[1][1] - Q[2][2] );

		C[11] = -Tau[10]*Q[0][1];

		C[12] = -Tau[10]*Q[1][2];

		C[13] = -Tau[10]*Q[0][2];

		C[14] = 0;

		



	    }


	
	    // Boundary node
	    
	    else {

		for( uint k = 0 ; k < q ; k++ ) {

		    C[k] = 0;

		}

	    }



	    

	}

	else {

	    cout << " [ERROR]  Surface tension term undefined" << endl << endl;

	    exit(1);

	}

    }

    

}



/** Model constant value */

const scalar liSurfaceTension::kappa() const {

    return _kappa;

}
