#include <pseudoPotEquation.H>

using namespace std;


/** Compute local density */

scalar pseudoPotEquation::localDensity( const uint& id ) {

    scalar r(0);

    const uint q = mesh.lmodel()->q();
    
    for( uint k = 0 ; k < q ; k++ )
    	r += _pdf[id][k];


    return r;    

}


/** Compute local velocity */
/** Density MUST be already updated */

const void pseudoPotEquation::localVelocity( std::vector<scalar>& v, const uint& id ) {
    

    // Local velocity
    
    vector<scalar> lv = {0,0,0};


    // Lattice constants

    const vector< vector<int> >& vel = mesh.lmodel()->lvel();

    const scalar q = mesh.lmodel()->q();


    // Compute first moment
    
    for( uint j = 0 ; j < 3 ; j++ ) {

    	for( uint k = 0 ; k < q ; k++ ) {

    	    lv[j] += vel[k][j] * _pdf[id][k];
		    
    	}
	    
    }


    // Add interaction force and divide by density

    vector<scalar> Ft = F.total(id);
    
    for( uint j = 0 ; j < 3 ; j++ )
    	lv[j] = ( lv[j]   +   0.5 * Ft[j]   ) / rho[id];
	


    // Copy to global array
    for( uint j = 0 ; j < 3 ; j++ ) 
    	v[j] = lv[j];
	

}









/** Constructor */

pseudoPotEquation::pseudoPotEquation( const string& name,
				      const latticeMesh& mesh_,
				      timeOptions& Time_,
				      pdfField& pdf_,
				      scalarField& rho_,
				      vectorField& U_,
				      scalarField& T_)
    
    : lbEquation(name, mesh_, Time_, pdf_),
      rho(rho_),
      U(U_),
      T(T_),
      F("properties/macroProperties", name, mesh_, Time_){    

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
