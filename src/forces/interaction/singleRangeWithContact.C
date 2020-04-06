#include <singleRangeWithContact.H>

using namespace std;


/** Constructor */

singleRangeWithContact::singleRangeWithContact( const string& dictName,
						const string& eqName,
						const latticeMesh& mesh,
						timeOptions& Time )

    : singleRangeIntForce(dictName, eqName, mesh, Time) {


    
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

singleRangeWithContact::~singleRangeWithContact() {}




/** Update force field. Use of virtual nodes if requested */

void singleRangeWithContact::update( scalarField& rho, scalarField& T ) {



    // Compute on bulk nodes first

    singleRangeIntForce::update( rho, T );

    
    
    // Reference to neighbour array

    const vector< vector<int> >& nb = _mesh.nbArray();


    // Lattice model properties

    vector< vector<int> > vel = _mesh.lmodel()->lvel();

    vector<uint> reverse = _mesh.lmodel()->reverse();

    scalar cs2 = _mesh.lmodel()->cs2();

    scalar q = _mesh.lmodel()->q();



    // Temporary

    const timeOptions& Time = rho.time();

    scalar cline[2] = {100000000, 0};

   

    // Move over points

    for( uint i = 0 ; i < _mesh.local() ; i++ ) {


	// Compute only if node is on boundary

	if( _mesh.isOnBoundary(i) ) {
	    


	    // Acum force

	    vector<scalar> F = {0, 0, 0};



	    // Compute aparent angle

	    scalar apangle( M_PI );

	    bool wallForce( false );

	    

	    if(      ( _withGeomContact )
	    	 &&  (  _contactAngle.find(i) != _contactAngle.end() )   ){ 
		
	    	// scalar gradRho[3] = {0,0,0};

	    	// rho.grad(gradRho, i);

	    	scalar gradRho[3] = { 0.5*(rho.at(nb[i][2]) - rho.at(nb[i][1])),
	    			      0.5*(rho.at(nb[i][4]) - rho.at(nb[i][3])),
	    			          // (rho.at(nb[i][6]) - rho.at(i))  };
	    			      0.5*( -3.0*rho.at(i) + 4.0*rho.at(nb[i][6]) - rho.at(nb[nb[i][6]][6]) )  };				      
				
	    	scalar gmag = sqrt( gradRho[0]*gradRho[0] + gradRho[1]*gradRho[1] );


	    	// Compute aparent angle

	    	if(gmag != 0) {
				
	    	    apangle = -gradRho[2]  /  gmag;

	    	    apangle = M_PI/2 - atan(apangle); 
			
	    	}


	    	if( apangle < _limitAngle )
	    	    wallForce = true;


		// if(wallForce)
		// {

		//     const vector<int>& point = _mesh.latticePoint(i);
		    
		//     scalar rad = sqrt( (point[0] - 30) * (point[0] - 30)
		// 		       + (point[1] - 30) * (point[1] - 30)
		// 		       + (point[2] - 0)  * (point[2] - 0) );


		//     if(rad > 8)
		// 	wallForce=false;



		// }		


	    }






	    if( wallForce ) {

		if(  _contactAngle.find(i) != _contactAngle.end() ){

		    const vector<int>& point = _mesh.latticePoint(i);
		    
		    scalar rad = sqrt( (point[0] - 30) * (point[0] - 30)
				       + (point[1] - 30) * (point[1] - 30)
				       + (point[2] - 0)  * (point[2] - 0) );

		    if(rad <= cline[0])
			cline[0] = rad;

		    if(rad >= cline[1])
			cline[1] = rad;

		}
		


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

				    // scalar gradRho[3] = {0,0,0};

				    // rho.grad(gradRho, first);

				    scalar gradRho[3] = { 0.5*(rho.at(nb[first][2]) - rho.at(nb[first][1])),
							  0.5*(rho.at(nb[first][4]) - rho.at(nb[first][3])),
							  0};	    
				
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


	    } // wallForce

		
	    

	}
	

    }





    if(_mesh.pid()==0){
    	if( Time.write(false) )
    	    cout << cline[0] << " " << cline[1] << endl;
    }
    

    // Sync across processors

    _force.sync();

}
