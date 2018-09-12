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






/** Default constructor */

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
