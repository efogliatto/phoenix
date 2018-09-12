#include <LiMRTEq.H>

using namespace std;


/** Default constructor */

LiMRTEq::LiMRTEq( const string& name,
		  const latticeMesh& mesh_,
		  timeOptions& Time_,
		  pdfField pdf_,
		  scalarField& rho_,
		  vectorField& U_,
		  scalarField& T_) : pseudoPotEquation(name,
						       mesh_,
						       Time_,
						       pdf_,
						       rho_,
						       U_,
						       T_) {


    for(uint i = 0 ; i < mesh.npoints() ; i++)
	cout << pdf_[i][0] << endl;
    
    // // Set equilibrium pdf values

    // LiMRTEq::setEquilibrium();

}
    

/** Default destructor */

LiMRTEq::~LiMRTEq() {}



/** Collision process */

const void LiMRTEq::collision() {}

    
/** Set pdf to equilibrium values */

const void LiMRTEq::setEquilibrium() {

    
    // Lattice model properties

    vector< vector<int> > vel = mesh.lmodel()->lvel();

    vector<scalar> omega = mesh.lmodel()->omega();
    
    scalar cs2 = mesh.lmodel()->cs2();
    

    // Compute equilibrium for all points
    
    for( uint i = 0 ; i < mesh.npoints() ; i++ ) {

    
    	for( uint k = 0 ; k < _q ; k++ ) {

    	    scalar alpha = 0,
    		beta = 0;


    	    // Dot product
	
    	    for( uint j = 0 ; j < 3 ; j++ ) {

    		alpha += vel[k][j] * U[i][j];

    		beta += U[i][j] * U[i][j];

    	    }

	    
    	    // pdf[i][k] = rho[i] * omega[k] * ( 1 + alpha/cs2   +   0.5 * alpha * alpha / (cs2*cs2)  -  0.5 * beta / cs2 );
	    // cout << pdf[i][0] << endl;	

    	}
	

    }
    

}
