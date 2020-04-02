#include <singleRangeMixedIntForce.H>

using namespace std;


/** Constructor */

singleRangeMixedIntForce::singleRangeMixedIntForce( const string& dictName,
						    const string& eqName,
						    const latticeMesh& mesh,
						    timeOptions& Time )

    : singleRangeIntForce(dictName, eqName, mesh, Time) {



    // Read beta constant
    
    dictionary dict( dictName );

    _beta = dict.lookUp<scalar>( eqName + "/Forces/Interaction/beta" );

    

}




/** Destructor */

singleRangeMixedIntForce::~singleRangeMixedIntForce() {}



/** Update force field */

void singleRangeMixedIntForce::update( scalarField& rho, scalarField& T ) {

    
    // Reference to neighbour array

    const vector< vector<int> >& nb = _mesh.nbArray();


    // Lattice model properties

    vector< vector<int> > vel = _mesh.lmodel()->lvel();

    vector<uint> reverse = _mesh.lmodel()->reverse();

    scalar cs2 = _mesh.lmodel()->cs2();

    scalar q = _mesh.lmodel()->q();
    

    // Move over points

    for( uint i = 0 ; i < _mesh.local() ; i++ ) {


	// Check if neighbour is on boundary

	bool isOnBnd(false);

	for( uint k = 0 ; k < q ; k++ ) {

	    if( nb[i][k] == -1 ) {

		isOnBnd = true;

		k = q;

	    }

	}


	


	// Compute only if node is not on boundary

	if( isOnBnd ) {


	    // No not compute flag

	    bool wallForce( true );

	    

	    if( _computeOnBnd ) {


		// Acum force

		vector<scalar> F1 = {0, 0, 0};

		vector<scalar> F2 = {0, 0, 0};



		// Compute aparent angle

		scalar apangle( M_PI );

		if( _withGeomContact ) {

		    scalar gradRho[3] = {0,0,0};

		    rho.grad(gradRho, i);
				
		    scalar gmag = sqrt( gradRho[0]*gradRho[0] + gradRho[1]*gradRho[1] );


		    // Compute aparent angle

		    if(gmag != 0) {
				
			apangle = -gradRho[2]  /  gmag;

			apangle = M_PI/2 - atan(apangle); 
			
		    }


		    if( apangle > (0.65 * M_PI) )
		    	wallForce = false;

		}
		


		// Move over velocities
		
		for( uint k = 1 ; k < q ; k++ ) {


		    // Neighbour index

		    int neighId = nb[i][ reverse[k] ];


		    // Virtual node

		    scalar _rho(0), _T(0);


		    if( neighId == -1 ) {


			// With geometric contact
		    
			if( _withGeomContact ) {



			    // Compute only if node index is on angles list

			    if( _contactAngle.find(i) != _contactAngle.end() ) {


				// Normal nodes

				int first = _mesh.vnode(i,k);

				int second = _mesh.vnode(i,k,false);


				

				// Compute apparent angle first
				
				// if( _hysteresis.at(i)[0] != _hysteresis.at(i)[1] ) {
			    					
				//     _contactAngle.at(i) = apangle;
				    
				// }



				/*********************************************/
				/*             Explicit scheme               */
				/*********************************************/
				
				
				// We need density gradients along the boundary

				scalar gradRho[3] = {0,0,0};

				rho.grad(gradRho, first);
				
				scalar gmag = sqrt( gradRho[0]*gradRho[0] + gradRho[1]*gradRho[1] );

				


				// Two point derivative
				
				_rho = rho.at(second) + 2*tan( M_PI/2 - _contactAngle.at(i) ) * gmag;

				_T = T.at( first );










				


				// /*********************************************/
				// /*             Implicit scheme               */
				// /*********************************************/
				
				
				// // We need density gradients along the boundary


				// // Initial gradient

				// scalar gradRho[3] = {0,0,0};

				// rho.grad(gradRho, first);												

				// scalar gmag = sqrt( gradRho[0]*gradRho[0] + gradRho[1]*gradRho[1] + gradRho[2]*gradRho[2] );



				// // Implicit calculation 

				// scalar vrho = rho.at(second) + gmag * cos(_contactAngle.at(first));

				// for(uint m = 0 ; m < 50 ; m++ ) {

				//     gradRho[2] = 0.5*( rho.at(nb[first][6]) - vrho );

				//     gmag = sqrt( gradRho[0]*gradRho[0] + gradRho[1]*gradRho[1] + gradRho[2]*gradRho[2] );

				//     vrho = rho.at(second) + gmag * cos(_contactAngle.at(first));

				// }

			       

				// _rho = vrho;

				// _T = T.at( first );
				


			    }

			    else {

				neighId = _mesh.vnode(i,k);

				_rho = rho.at( neighId );
			    
				_T = T.at( neighId );

			    }

			}

			
			else {

			    neighId = _mesh.vnode(i,k);

			    _rho = rho.at( neighId );
			    
			    _T = T.at( neighId );			    

			}
			

		    }

		    else {

			_rho = rho.at( neighId );
			    
			_T = T.at( neighId );

		    }


		    
		    

		    // Complete force

		    scalar alpha_1 = _weights[k] * potential( _rho, _T, cs2 );

		    scalar alpha_2 = alpha_1 * alpha_1 / _weights[k];

	    
		    for( uint j = 0 ; j < 3 ; j++ ) {

			F1[j] +=  alpha_1 * (scalar)vel[k][j] ;

			F2[j] +=  alpha_2 * (scalar)vel[k][j] ;

		    }

		    
		    

		}		



		// Extra constant
		
		scalar kappa = -_G * potential( rho[i], T[i], cs2 ) * _beta;
    
		for( uint j = 0 ; j < 3 ; j++) {
	
		    _force[i][j] =  F1[j] * kappa
 			         +   0.5 * ( 1 - _beta) * F2[j] ;
	
		}







		// Correct force 
		
		if( wallForce == false ) {

		    for( uint j = 0 ; j < 3 ; j++ )
			_force[i][j] = 0;


		}

	    }




	    // Do not compute on boundary if not enabled
	    
	    else {
	    
		for( uint j = 0 ; j < 3 ; j++ )
			_force[i][j] = 0;
		
	    }
	    

	}


	

	// Bulk fluid node
	
	else {

	    
	    vector<scalar> F1 = {0, 0, 0};

	    vector<scalar> F2 = {0, 0, 0};	    

	    for( uint k = 1 ; k < q ; k++ ) {
       
		int neighId = nb[i][ reverse[k] ];

		scalar _rho = rho[neighId];

		scalar _T = T[neighId];
	    
		scalar alpha_1 = _weights[k] * potential( _rho, _T, cs2 );

		scalar alpha_2 = alpha_1 * alpha_1 / _weights[k];

	    
		for( uint j = 0 ; j < 3 ; j++ ) {

		    F1[j] +=  alpha_1 * (scalar)vel[k][j] ;

		    F2[j] +=  alpha_2 * (scalar)vel[k][j] ;

		}
    

	    }

		

	    // Extra constant
		
	    scalar kappa = -_G * potential( rho[i], T[i], cs2 ) * _beta;
    
	    for( uint j = 0 ; j < 3 ; j++) {
	
		_force[i][j] =  F1[j] * kappa
		              + 0.5 * ( 1 - _beta) * F2[j] ;
	
	    }	

	    

	}
	

    }




    // Sync across processors

    _force.sync();

}






// /** Update force field */

// void singleRangeMixedIntForce::update( scalarField& rho, scalarField& T ) {

    
//     // Reference to neighbour array

//     const vector< vector<int> >& nb = _mesh.nbArray();


//     // Lattice model properties

//     vector< vector<int> > vel = _mesh.lmodel()->lvel();

//     vector<uint> reverse = _mesh.lmodel()->reverse();

//     scalar cs2 = _mesh.lmodel()->cs2();

//     scalar q = _mesh.lmodel()->q();
    

//     // Move over points

//     for( uint i = 0 ; i < _mesh.local() ; i++ ) {


// 	// Check if neighbour is on boundary

// 	bool isOnBnd(false);

// 	for( uint k = 0 ; k < q ; k++ ) {

// 	    if( nb[i][k] == -1 ) {

// 		isOnBnd = true;

// 		k = q;

// 	    }

// 	}



// 	// Compute only if node is not on bounday

// 	if( isOnBnd ) {

// 	    for( uint j = 0 ; j < 3 ; j++ )
// 		_force[i][j] = 0;

// 	}


// 	else {

	    
// 	    vector<scalar> F1 = {0, 0, 0};

// 	    vector<scalar> F2 = {0, 0, 0};	    

// 	    for( uint k = 1 ; k < q ; k++ ) {
       
// 		int neighId = nb[i][ reverse[k] ];

// 		scalar _rho = rho[neighId];

// 		scalar _T = T[neighId];
	    
// 		scalar alpha_1 = _weights[k] * potential( _rho, _T, cs2 );

// 		scalar alpha_2 = alpha_1 * alpha_1 / _weights[k];

	    
// 		for( uint j = 0 ; j < 3 ; j++ ) {

// 		    F1[j] +=  alpha_1 * (scalar)vel[k][j] ;

// 		    F2[j] +=  alpha_2 * (scalar)vel[k][j] ;

// 		}
    

// 	    }

		

// 	    // Extra constant
		
// 	    scalar kappa = -_G * potential( rho[i], T[i], cs2 ) * _beta;
    
// 	    for( uint j = 0 ; j < 3 ; j++) {
	
// 		_force[i][j] =  F1[j] * kappa
// 		              + 0.5 * ( 1 - _beta) * F2[j] ;
	
// 	    }	

	    

// 	}
	

//     }




//     // Sync across processors

//     _force.sync();

// }
