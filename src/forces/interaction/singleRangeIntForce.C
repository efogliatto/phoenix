#include <singleRangeIntForce.H>

using namespace std;


/** Constructor */

singleRangeIntForce::singleRangeIntForce( const string& dictName,
					  const string& eqName,
					  const latticeMesh& mesh,
					  timeOptions& Time )

    : interactionForce(dictName, eqName, mesh, Time) {



    // Update pseudo pot weights

    string model = mesh.lmodel()->name();

    if( model == "D2Q9" ) {

    	_weights.push_back( 0 );
    	_weights.push_back( 1.0/3 );
    	_weights.push_back( 1.0/3 );
    	_weights.push_back( 1.0/3 );
    	_weights.push_back( 1.0/3 );
    	_weights.push_back( 1.0/12 );
    	_weights.push_back( 1.0/12 );
    	_weights.push_back( 1.0/12 );
    	_weights.push_back( 1.0/12 );

    }

    else {

	if( model == "D3Q15" ) {

	    _weights.push_back( 0 );

	    for( uint i = 1 ; i <= 6 ; i++ )
		_weights.push_back( 1.0 / 3.0 );

	    for( uint i = 7 ; i < 15 ; i++ )
		_weights.push_back( 1.0 / 24.0 );
	    


	}

	else {

	    cout << " [ERROR]  Potential weights not yet implemented for " << model << endl;

	}

    }

}




/** Destructor */

singleRangeIntForce::~singleRangeIntForce() {}





/** Update force field. Use of virtual nodes if requested */

void singleRangeIntForce::update( scalarField& rho, scalarField& T ) {

    
    
    // Reference to neighbour array

    const vector< vector<int> >& nb = _mesh.nbArray();


    // Lattice model properties

    vector< vector<int> > vel = _mesh.lmodel()->lvel();

    vector<uint> reverse = _mesh.lmodel()->reverse();

    scalar cs2 = _mesh.lmodel()->cs2();

    scalar q = _mesh.lmodel()->q();


    // Apparent angle
    
    // scalar apangle = apparentContactAngle(rho, "Y0");
    // scalar apangle(0);

    

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

	    if( _computeOnBnd ) {


		// Acum force

		vector<scalar> F = {0, 0, 0};
		


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


				

				// // Compute apparent angle first
				
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

		    scalar alpha = _weights[k] * potential( _rho, _T, cs2 );
		    
	    
		    for( uint j = 0 ; j < 3 ; j++ ) {

		    	F[j] +=  alpha * (scalar)vel[k][j] ;

		    }		    
		    

		}		



		// Extra constant
		
		scalar beta = -_G * potential( rho.at(i), T.at(i), cs2 );
    
		for( uint j = 0 ; j < 3 ; j++) {
	
		    _force[i][j] =  F[j] * beta;
	
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

	    
	    vector<scalar> F = {0, 0, 0};	    

	    for( uint k = 1 ; k < q ; k++ ) {
       
		int neighId = nb[i][ reverse[k] ];

		scalar _rho = rho[neighId];

		scalar _T = T[neighId];
	    
		scalar alpha = _weights[k] * potential( _rho, _T, cs2 );

	    
		for( uint j = 0 ; j < 3 ; j++ ) {

		    F[j] +=  alpha * (scalar)vel[k][j] ;

		}
    

	    }

		

	    // Extra constant
		
	    scalar beta = -_G * potential( rho[i], T[i], cs2 );
    
	    for( uint j = 0 ; j < 3 ; j++) {
	
		_force[i][j] =  F[j] * beta;
	
	    }	

	    

	}
	

    }




    // Sync across processors

    _force.sync();

}
