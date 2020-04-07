#include <liSurfaceTensionContact.H>

#include <dictionary.H>

using namespace std;


/** Constructor */

liSurfaceTensionContact::liSurfaceTensionContact( const string& dictName,
						  const string& eqName,
						  const latticeMesh& mesh )
    
    : liSurfaceTension(dictName, eqName, mesh) {


    
    // // Update model constant
    
    // dictionary dict(dictName);

    // _kappa = dict.lookUp<scalar>( eqName + "/SurfaceTensionContact/kappa" );
    

}




/** Destructor */

liSurfaceTensionContact::~liSurfaceTensionContact() {}






/** Force at specific node */

const void liSurfaceTensionContact::ST( const uint& i, const scalarField& rho, const scalarField& T, vector<scalar>& C, interactionForce* _fi, const vector<scalar>& Tau ) const {


    // Base version for bulk node

    liSurfaceTension::ST(i, rho, T, C, _fi, Tau);


    // if( _mesh.isOnBoundary(i) ) {

	
    // 	// Lattice constants
    
    // 	const uint q = _mesh.lmodel()->q();

    // 	const scalar cs2 = _mesh.lmodel()->cs2();

    // 	const vector< vector<int> > vel = _mesh.lmodel()->lvel();

    // 	const vector<uint>& reverse = _mesh.lmodel()->reverse();

    // 	const vector< vector<int> >& nb = _mesh.nbArray();



    // 	if( q == 15 ) {


    // 	    // Compute aparent angle

    // 	    scalar apangle( M_PI );

    // 	    bool wallForce( false );
    

    // 	    if(  _contactAngle.find(i) != _contactAngle.end()  ){ 
		
    // 	    	scalar gradRho[3] = {0,0,0};

    // 		rho.cartesianGradient(gradRho, i);		
				
    // 	    	scalar gmag = sqrt( gradRho[0]*gradRho[0] + gradRho[1]*gradRho[1] );


    // 	    	// Compute aparent angle

    // 	    	if(gmag != 0) {
				
    // 	    	    apangle = -gradRho[2]  /  gmag;

    // 	    	    apangle = M_PI/2 - atan(apangle); 
			
    // 	    	}


    // 	    	if( apangle < _limitAngle )
    // 	    	    wallForce = true;


    // 	    }


    // 	    if( wallForce ) {

	    
	    
    // 		// Q tensor

    // 		scalar Q[3][3] = { {0,0,0}, {0,0,0}, {0,0,0} };

    // 		const scalar phi = _fi->potential(rho.at(i), T.at(i), cs2);

    // 		const scalar omega[15] = { 0, 1.0/3, 1.0/3, 1.0/3, 1.0/3, 1.0/3, 1.0/3, 1.0/24, 1.0/24, 1.0/24, 1.0/24, 1.0/24, 1.0/24, 1.0/24, 1.0/24  };
	


    // 		// Move over velocities

    // 		for( uint k = 1 ; k < q ; k++ ) {

    // 		    int neighId = nb[i][ reverse[k] ];

    // 		    scalar alpha(0);


    // 		    // Existing neighbours
		
    // 		    if (neighId != -1) {

    // 			alpha = omega[k] * ( _fi->potential( rho.at(neighId), T.at(neighId), cs2 ) - phi );

    // 		    }


    // 		    // Ghost neighbour
		
    // 		    else {


    // 			// Normal nodes

    // 			int first = _mesh.vnode(i,k);
			
    // 			int second = _mesh.vnode(i,k,false);


    // 			// Compute only if node index is on angles list

    // 			if( _contactAngle.find(first) != _contactAngle.end() ) {
			    

    // 			    /*********************************************/
    // 			    /*             Explicit scheme               */
    // 			    /*********************************************/
				
				
    // 			    // We need density gradients along the boundary

    // 			    scalar gradRho[3] = {0,0,0};

    // 			    rho.cartesianGradient(gradRho, first);			    
				
    // 			    scalar gmag = sqrt( gradRho[0]*gradRho[0] + gradRho[1]*gradRho[1] );
				


    // 			    // Two point derivative
				
    // 			    scalar _rho = rho.at(second) + 2*tan( M_PI/2 - _contactAngle.at(i) ) * gmag;

    // 			    scalar _T = T.at( first );

    // 			    alpha = omega[k] * ( _fi->potential( _rho, _T, cs2 ) - phi );
			       

    // 			}

    // 			else {

    // 			    scalar _rho = rho.at( first );
			    
    // 			    scalar _T = T.at( first );

    // 			    alpha = omega[k] * ( _fi->potential( _rho, _T, cs2 ) - phi );

    // 			}			

			

    // 		    }
		


		


    // 		    for( uint ii = 0 ; ii < 3 ; ii++ ) {

    // 			for( uint jj = 0 ; jj < 3 ; jj++ ) {

    // 			    Q[ii][jj] += alpha * vel[k][ii] * vel[k][jj];

    // 			}

    // 		    }
		    

    // 		}


    // 		for( uint ii = 0 ; ii < 3 ; ii++ ) {

    // 		    for( uint jj = 0 ; jj < 3 ; jj++ ) {

    // 			Q[ii][jj] = -Q[ii][jj] * 0.5 * phi * _kappa;

    // 		    }

    // 		}




    // 		// Complete C array

    // 		C[0] = 0;

    // 		C[1] = (4.0/5.0) * Tau[1] * (Q[0][0] + Q[1][1] + Q[2][2]);

    // 		C[2] = 0;

    // 		C[3] = 0;

    // 		C[4] = 0;

    // 		C[5] = 0;

    // 		C[6] = 0;

    // 		C[7] = 0;

    // 		C[8] = 0;

    // 		C[9] = -Tau[10]*( 2.0*Q[0][0] - Q[1][1] - Q[2][2] );

    // 		C[10] = -Tau[10]*( Q[1][1] - Q[2][2] );

    // 		C[11] = -Tau[10]*Q[0][1];

    // 		C[12] = -Tau[10]*Q[1][2];

    // 		C[13] = -Tau[10]*Q[0][2];

    // 		C[14] = 0;







    // 	    } // wallForce
	    

    // 	}
	

    // }


    

}





/** Read geometric contact angle properties */

const void liSurfaceTensionContact::readGeometricContact( const string& eqname ) {

    
    // // Read contact angles for each boundary

    // dictionary dict("properties/macroProperties");

    // string contact = dict.lookUpOrDefault<string>( eqname + "/Forces/Interaction/ContactAngle/type", "none");
    

    // if( contact == "geometric" ) {



    // 	// Limit apparent angle

    // 	_limitAngle = dict.lookUpOrDefault<scalar>( eqname + "/Forces/Interaction/ContactAngle/limitAngle", 180 );

    // 	_limitAngle = _limitAngle * M_PI / 180.0;	

	

    // 	// Read angles for each boundary

    // 	const map< string, vector<uint> >& boundary = _mesh.boundaries();

    // 	map<string, scalar> angle;

    // 	map<string, scalar> hmin;

    // 	map<string, scalar> hmax;
	

    // 	for( const auto& bd : boundary) {

    // 	    scalar ang = dict.lookUpOrDefault<scalar>( eqname + "/Forces/Interaction/ContactAngle/" + bd.first + "/Theta/", -1);

    // 	    if(ang >= 0) {
		
    // 		angle[bd.first] = ang * M_PI / 180.0;


    // 		hmin[bd.first] = dict.lookUpOrDefault<scalar>( eqname + "/Forces/Interaction/ContactAngle/" + bd.first + "/Hysteresis/min", 0);

    // 		hmax[bd.first] = dict.lookUpOrDefault<scalar>( eqname + "/Forces/Interaction/ContactAngle/" + bd.first + "/Hysteresis/max", 0);		
		

    // 	    }

    // 	}



    // 	// Create map for each boundary node

    // 	for( const auto& bd : angle) {

    // 	    for( const auto& id : boundary.at(bd.first) ) {

    // 		_contactAngle[id] = bd.second;

    // 		_hysteresis[id] = {hmin[bd.first], hmax[bd.first]};

    // 	    }

    // 	}	
	

    // }

}
