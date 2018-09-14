#include <LiMRTEq.H>

using namespace std;


/** Default constructor */

LiMRTEq::LiMRTEq( const string& name,
		  const latticeMesh& mesh_,
		  timeOptions& Time_,
		  pdfField& pdf_,
		  scalarField& rho_,
		  vectorField& U_,
		  scalarField& T_)

    // Base class construction
    
    : pseudoPotEquation(name,
			mesh_,
			Time_,
			pdf_,
			rho_,
			U_,
			T_) {


    
    // Set equilibrium pdf values

    LiMRTEq::setEquilibrium();

}
    

/** Default destructor */

LiMRTEq::~LiMRTEq() {}




/** Collision process */

const void LiMRTEq::collision() {

    
    // Lattice constants

    scalarMatrix M = mesh.lmodel()->MRTMatrix();

    scalarMatrix invM = mesh.lmodel()->MRTInvMatrix();    

    const uint nodes = mesh.npoints();


    
    // Partial distributions
    
    vector<scalar> m     = {0,0,0,0,0,0,0,0,0};   // m:  momentum space
    
    vector<scalar> m_eq  = {0,0,0,0,0,0,0,0,0};   // meq: equilibrium in momentum space

    vector<scalar> S     = {0,0,0,0,0,0,0,0,0};   // MRT force

    vector<scalar> C     = {0,0,0,0,0,0,0,0,0};   // Surface tension term


    

    // Move over all points

    for( uint id = 0 ; id < nodes ; id++ ) {


	// Velocity magnitude

	scalar umag (0);
	
	for( uint j = 0 ; j < 3 ; j++ )	
	    umag += U[id][j] * U[id][j];


	
	// Compute equilibrium in momentum space
	
	m_eq[0] = rho[id];
	m_eq[1] = rho[id] * (-2 + 3*umag);
	m_eq[2] = rho[id] * (1 - 3*umag);
	m_eq[3] = rho[id] *   U[id][0];
	m_eq[4] = rho[id] * (-U[id][0]);
	m_eq[5] = rho[id] *   U[id][1];
	m_eq[6] = rho[id] * (-U[id][1]);
	m_eq[7] = rho[id] * (U[id][0]*U[id][0] - U[id][1]*U[id][1]);
	m_eq[8] = rho[id] * U[id][0] * U[id][1];
	

	// Distribution in momentum space

	// M.matDotVec(_pdf[id], m_eq);
	

    }

    

}




/** Set pdf to equilibrium values */

const void LiMRTEq::setEquilibrium() {

    
    // Lattice model properties

    vector< vector<int> > vel = mesh.lmodel()->lvel();

    vector<scalar> omega = mesh.lmodel()->omega();
    
    scalar cs2 = mesh.lmodel()->cs2();

    const uint np = mesh.npoints();

    const uint q = mesh.lmodel()->q();
    

    // Compute equilibrium for all points
    
    for( uint i = 0 ; i < np ; i++ ) {

    
    	for( uint k = 0 ; k < q ; k++ ) {

    	    scalar alpha = 0,
    	    	beta = 0;


    	    // Dot product
	
    	    for( uint j = 0 ; j < 3 ; j++ ) {

    	    	alpha += vel[k][j] * U[i][j];

    	    	beta += U[i][j] * U[i][j];

    	    }

	    
    	    _pdf.set(i, k, rho[i] * omega[k] * ( 1 + alpha/cs2   +   0.5 * alpha * alpha / (cs2*cs2)  -  0.5 * beta / cs2 ) );

    	}
	

    }
    

}
