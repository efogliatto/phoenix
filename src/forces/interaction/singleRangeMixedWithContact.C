#include <singleRangeMixedWithContact.H>

using namespace std;



/** Constructor */

singleRangeMixedWithContact::singleRangeMixedWithContact( const string& dictName,
							  const string& eqName,
							  const latticeMesh& mesh,
							  timeOptions& Time )

    : singleRangeMixedIntForce(dictName, eqName, mesh, Time) {


    
    // Read beta constant
    
    dictionary dict( dictName );

    string scheme = dict.lookUpOrDefault<string>( eqName + "/Forces/Interaction/ContactAngle/scheme", "explicit" );



    // geometrical schemes

    map<string, contactType> _cMapType;

    _cMapType["explicit"] = contactType::exp;

    _cMapType["implicit"] = contactType::imp;


    
    // Assign scheme
    
    if( _cMapType.find(scheme) != _cMapType.end() ) {
	

	switch( _cMapType[scheme] ) {
	

	case contactType::exp:

	    _scheme = contactType::exp;

	    break;


	case contactType::imp:

	    _scheme = contactType::imp;

	    break;	    

	}
	

    }

    else {

	cout << endl << " [ERROR]  Geometrical scheme " << scheme << " not available" << endl << endl;

	exit(1);

    }






    // Limit apparent angle

    _limitAngle = dict.lookUpOrDefault<scalar>( eqName + "/Forces/Interaction/ContactAngle/limitAngle", 180 );

    _limitAngle = _limitAngle * M_PI / 180.0;
    

}




/** Destructor */

singleRangeMixedWithContact::~singleRangeMixedWithContact() {}




/** Update force field. Use of virtual nodes if requested */

void singleRangeMixedWithContact::update( scalarField& rho, scalarField& T ) {



    // Compute on bulk nodes first

    singleRangeMixedIntForce::update( rho, T );

    
    
    // Reference to neighbour array

    const vector< vector<int> >& nb = _mesh.nbArray();


    // Lattice model properties

    vector< vector<int> > vel = _mesh.lmodel()->lvel();

    vector<uint> reverse = _mesh.lmodel()->reverse();

    scalar cs2 = _mesh.lmodel()->cs2();

    scalar q = _mesh.lmodel()->q();



   

    // Move over points

    for( uint i = 0 ; i < _mesh.local() ; i++ ) {


	// Compute only if node is not on boundary

	if( _mesh.isOnBoundary(i) ) {
	    


	    // Acum force

	    vector<scalar> F1 = {0, 0, 0};

	    vector<scalar> F2 = {0, 0, 0};




	    // Compute aparent angle

	    scalar apangle( M_PI );

	    bool wallForce( true );

	    

	    if( _withGeomContact ) {
		
		scalar gradRho[3] = {0,0,0};

		rho.grad(gradRho, i);
				
		scalar gmag = sqrt( gradRho[0]*gradRho[0] + gradRho[1]*gradRho[1] );


		// Compute aparent angle

		if(gmag != 0) {
				
		    apangle = -gradRho[2]  /  gmag;

		    apangle = M_PI/2 - atan(apangle); 
			
		}


		if( apangle > _limitAngle )
		    wallForce = false;


	    }



	    if( wallForce ) {
		


		// Move over velocities
		
		for( uint k = 1 ; k < q ; k++ ) {


		    // Neighbour index

		    int neighId = nb[i][ reverse[k] ];


		    // Virtual node

		    scalar _rho(0), _T(0);


		    if( neighId == -1 ) {



			// Normal nodes

			int first = _mesh.vnode(i,k);

			int second = _mesh.vnode(i,k,false);
		    


			// With geometric contact
		    
			if( _withGeomContact ) {



			    // Compute only if node index is on angles list

			    if( _contactAngle.find(first) != _contactAngle.end() ) {
			    

				
				if( _scheme == contactType::exp ) {



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



				}


			    
			    
				if( _scheme == contactType::imp ) {



				    /*********************************************/
				    /*             Implicit scheme               */
				    /*********************************************/
				
				
				    // We need density gradients along the boundary


				    // Initial gradient

				    scalar gradRho[3] = {0,0,0};

				    rho.grad(gradRho, first);												

				    scalar gmag = sqrt( gradRho[0]*gradRho[0] + gradRho[1]*gradRho[1] + gradRho[2]*gradRho[2] );



				    // Implicit calculation 

				    scalar vrho = rho.at(second) + gmag * cos(_contactAngle.at(first));

				    for(uint m = 0 ; m < 50 ; m++ ) {

					gradRho[2] = 0.5*( rho.at(nb[first][6]) - vrho );

					gmag = sqrt( gradRho[0]*gradRho[0] + gradRho[1]*gradRho[1] + gradRho[2]*gradRho[2] );

					vrho = rho.at(second) + gmag * cos(_contactAngle.at(first));

				    }

			       

				    _rho = vrho;

				    _T = T.at( first );


				}
				


			    }

			    else {

				neighId = _mesh.vnode(i,k);

				_rho = rho.at( neighId );
			    
				_T = T.at( neighId );

			    }

			}





			// No geometric model
		    
			else {

			    _rho = 2.0*rho.at( first ) - rho.at( second );
			    
			    _T = 2.0*T.at( first ) - T.at( second );

			}
		    
			

		    }




		    // Neighbour exists
		
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


	    } // wallForce

		
	    

	}
	

    }




    // Sync across processors

    _force.sync();

}
