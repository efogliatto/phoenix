#include <pseudoPotEquation.H>

using namespace std;


/** Compute local density */

scalar pseudoPotEquation::localDensity( const uint& id ) {

    scalar r(0);
	
    for( uint k = 0 ; k < _q ; k++ )
    	r += pdf[id][k];


    return r;    

}


/** Compute local velocity */
/** Density MUST be already updated */

const void pseudoPotEquation::localVelocity( std::vector<scalar>& v, const uint& id ) {
    

    // Local velocity
    
    scalar lv[3] = {0,0,0};


    // Lattice velocities

    const vector< vector<int> > vel = mesh.lmodel()->lvel();


    // Compute first moment
    
    for( uint j = 0 ; j < 3 ; j++ ) {

    	for( uint k = 0 ; k < _q ; k++ ) {

    	    lv[j] += vel[k][j] * pdf[id][k];
		    
    	}
	    
    }


    // Add interaction force and divide by density

    for( uint j = 0 ; j < 3 ; j++ )
    	lv[j] = ( lv[j]   +   0.5 * Ft[id][j]   ) / rho[id];
	


    // Copy to global array
    for( uint j = 0 ; j < 3 ; j++ ) 
    	v[j] = lv[j];
	

}










/** Interaction force field */

const void pseudoPotEquation::updateIntForce() {


    // Reference to neighbour array

    const vector< vector<int> > nb = mesh.nbArray();


    // Lattice model properties

    vector< vector<int> > vel = mesh.lmodel()->lvel();

    vector<uint> reverse = mesh.lmodel()->reverse();

    scalar cs2 = mesh.lmodel()->cs2();

    
    
    

    for( uint id = 0 ; id < mesh.local() ; id++ ) {
    
	vector<scalar> F = {0.0, 0.0, 0.0};

	uint noneigh(0);
    

	// Move over neighbours and check for boundary
	for( uint k = 1 ; k < _q ; k++ ) {

		if( nb[id][k] == -1 ) {
	    
		    noneigh++;
	    
		}

	}


	// Do not use unexisting neighbour
	if( noneigh == 0 ) {

	
		for( uint k = 1 ; k < _q ; k++ ) {
       

		    int neighId = nb[id][ reverse[k] ];

		    scalar nbrho = rho[neighId];

		    scalar nbt = T[neighId];
	    
		    scalar alpha = weights[k] * eos->potential( eos->p_eos(nbrho, nbt), cs2, nbrho, nbt );

	    
		    for( uint i = 0 ; i < 3 ; i++ ) {

		    	F[i] +=  alpha * (scalar)vel[k][i] ;

		    }
    

		}

	

	

		// Extra constant
		
		scalar beta = -eos->G() * eos->potential( eos->p_eos(rho[id], T[id]), cs2, rho[id], T[id] );
    
		for( uint i = 0 ; i < 3 ; i++) {
	
		    Fi[id][i] =  F[i] * beta;
	
		}	

		
	}

    }

}


    
/** Total force field */

const void pseudoPotEquation::updateTotalForce() {

    for( uint i = 0 ; i < mesh.local() ; i++ ) {
    
	for( uint j = 0 ; j < 3 ; j++) {
	
	    Ft[i][j] = Fi[i][j] +  (rho[i] - rho_0) * gravity[j];
	
	}

    }


}






/** Constructor */

pseudoPotEquation::pseudoPotEquation( const string& name,
				      const latticeMesh& mesh_,
				      timeOptions& Time_,
				      pdfField pdf_,
				      scalarField& rho_,
				      vectorField& U_,
				      scalarField& T_) : lbEquation(name, mesh_, Time_, pdf_),
							 rho(rho_),
							 U(U_),
							 T(T_),
							 Fi(mesh,Time,"Fi",IO::NO_READ,IO::NO_WRITE),
							 Ft(mesh,Time,"Ft",IO::NO_READ,IO::NO_WRITE)  {


    // Create EOS

    EOSCreator creator;

    eos = creator.create("properties/macroProperties");

    

    // Read gravity field and ext fixed force

    dictionary dict("properties/macroProperties");

    gravity  = dict.lookUp< vector<scalar> >("Navier-Stokes/gravity");

    extForce = dict.lookUp< vector<scalar> >("Navier-Stokes/extForce");


    
    // Create pp weights array

    if(   ( mesh.lmodel()->d() == 2 )  &&  ( mesh.lmodel()->q() == 9 )   ) {

	weights.push_back( 0 );
	weights.push_back( 1.0/3 );
	weights.push_back( 1.0/3 );
	weights.push_back( 1.0/3 );
	weights.push_back( 1.0/3 );
	weights.push_back( 1.0/12 );
	weights.push_back( 1.0/12 );
	weights.push_back( 1.0/12 );
	weights.push_back( 1.0/12 );	

    }

    else {

	cout << " [ERROR]  Pseudopotential model not yet implemented" << endl << endl;

    }
    

}


/** Default destructor */

pseudoPotEquation::~pseudoPotEquation() {}



/** Collision process */

const void pseudoPotEquation::collision() {}



/** Update macroscopic density */

const void pseudoPotEquation::updateMacroDensity() {

    for( uint i = 0 ; i < mesh.npoints() ; i++ )
	rho[i] = localDensity(i);

}

    
/** Update macroscopic velocity */

const void pseudoPotEquation::updateMacroVelocity() {

    for( uint i = 0 ; i < mesh.npoints() ; i++ )
	localVelocity( U[i], i );

}
