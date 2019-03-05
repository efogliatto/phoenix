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

const void liSurfaceTension::ST( const uint& i, const scalarField& rho, const scalarField& T, vector<scalar>& C, interactionForce* _fi, const vector<scalar>& Tau ) const {


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


	// Move over neighbours and check for boundary

	uint noneigh(0);

	for( uint k = 1 ; k < q ; k++ ) {

	    if( nb[i][k] == -1 ) {
	    
		noneigh++;
	    
	    }

	}


	// Do not use unexisting neighbour
	if( noneigh == 0 ) {


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












	

	// scalar psi = potential( mesh, mfields->rho[id], mfields->T[id]);

	// scalar alpha;

	// int neighId;



	    // // Move over neighbours and check for boundary

	    // for( k = 1 ; k < mesh->lattice.Q ; k++ ) {

	    // 	if( mesh->mesh.nb[id][k] == -1 ) {
	    
	    // 	    noneigh++;
	    
	    // 	}

	    // }	    


	    // // Do not use unexisting neighbour
	    // if( noneigh == 0 ) {

		
		// Q = matrixDoubleAlloc(2, 2, 0);


		// // Move over velocities

		// for( k = 1 ; k < 9 ; k++ ) {

		//     neighId = mesh->mesh.nb[id][ mesh->lattice.reverse[k] ];
		
		//     alpha = mesh->lattice.weights[k] * ( potential( mesh, mfields->rho[neighId], mfields->T[neighId]) - psi );


		//     for( i = 0 ; i < 2 ; i++ ) {

		// 	for( j = 0 ; j < 2 ; j++ ) {

		// 	    Q[i][j] += alpha * mesh->lattice.vel[k][i] * mesh->lattice.vel[k][j];

		// 	}

		//     }
		    

		// }


		// for( i = 0 ; i < 2 ; i++ ) {

		//     for( j = 0 ; j < 2 ; j++ ) {

		// 	Q[i][j] = Q[i][j] * 0.5 * psi * mesh->EOS.G * field->lbparam.liMRT.kappa_st;

		//     }

		// }



		// // Complete C array

		// C[0] = 0;

		// C[1] = 1.5 * field->lbparam.liMRT.Lambda[1] * (Q[0][0] + Q[1][1]);

		// C[2] = -1.5 * field->lbparam.liMRT.Lambda[2] * (Q[0][0] + Q[1][1]);

		// C[3] = 0;

		// C[4] = 0;

		// C[5] = 0;

		// C[6] = 0;

		// C[7] = -field->lbparam.liMRT.Lambda[8] * (Q[0][0] - Q[1][1]);

		// C[8] = -field->lbparam.liMRT.Lambda[8] * Q[0][1];






	    // }


	    // // Boundary node
	    
	    // else {

	    // 	for( k = 0 ; k < 9 ; k++ ) {

	    // 	    C[k] = 0;

	    // 	}

	    // }

	


	

    }

    else {

	if( q == 15 ) {



	}

	else {

	    cout << " [ERROR]  Surface tension term undefined" << endl << endl;

	    exit(1);

	}

    }
    
    
    // if( closestNodes.find(i) != closestNodes.end() ) {

	
    // 	// Reference to neighbour array

    // 	const vector< vector<int> >& nb = _mesh.nbArray();


    // 	// Lattice model properties

    // 	const vector< vector<int> >& vel = _mesh.lmodel()->lvel();

    // 	const vector<uint>& reverse = _mesh.lmodel()->reverse();

    // 	const scalar cs2 = _mesh.lmodel()->cs2();

    // 	const uint q = _mesh.lmodel()->q();


    // 	vector<scalar> F = {0, 0, 0};	    

    // 	for( uint k = 1 ; k < q ; k++ ) {
       
    // 	    int neighId = nb[i][ reverse[k] ];

    // 	    if( _mesh.isOnBoundary(neighId) ) {
	    
    // 		for( uint j = 0 ; j < 3 ; j++ ) {

    // 		    F[j] +=  _weights[k] * (scalar)vel[k][j] ;

    // 		}

    // 	    }
    

    // 	}

		

    // 	// Extra constant
		
    // 	scalar beta = potential( rho.at(i), T.at(i), cs2 );

    // 	beta = beta * closestNodes.at(i);	    
    
    // 	for( uint j = 0 ; j < 3 ; j++)
    // 	    f[j] =  F[j] * beta;
	
	

    // }

}
