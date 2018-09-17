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



    // Read sigma constant

    dictionary dict("properties/macroProperties");

    _sigma = dict.lookUp<scalar>( name + "/LBModel/sigma" );
    
    
    
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

    const scalar cs2 = mesh.lmodel()->cs2();

    const scalar q = mesh.lmodel()->q();

    scalar Ft[3];

    
    // Partial distributions
    
    vector<scalar> m(9);     // m:  momentum space
    
    vector<scalar> m_eq(9);  // meq: equilibrium in momentum space

    vector<scalar> S(9);     // MRT force

    // vector<scalar> C(9);     // Surface tension term


    

    // Move over all points

    for( uint id = 0 ; id < nodes ; id++ ) {


	// Local velocity and density

	const scalar r = rho[id];

	const scalar u[3] = { U.at(id,0), U.at(id,1), U.at(id,2) };
	

	// Velocity magnitude

	scalar umag (0);
	
	for( uint j = 0 ; j < 3 ; j++ )	
	    umag += u[j] * u[j];


	
	// Compute equilibrium in moment space
	
	m_eq[0] = r;
	m_eq[1] = r * (-2 + 3*umag);
	m_eq[2] = r * (1 - 3*umag);
	m_eq[3] = r *  u[0];
	m_eq[4] = r * (-u[0]);
	m_eq[5] = r *   u[1];
	m_eq[6] = r * (-u[1]);
	m_eq[7] = r * (u[0]*u[0] - u[1]*u[1]);
	m_eq[8] = r * u[0] * u[1];

	

	// Distribution in moment space

	M.matDotVec( _pdf[id], m );


	
	// Source term

	scalar psi = F.potential( r, T.at(id), cs2 );

	F.total(Ft, id);

	const scalar Fi[3] = { F.interaction(id,0), F.interaction(id,1), F.interaction(id,2) };


	S[0] = 0;
	S[1] =  6 * (u[0]*Ft[0] + u[1]*Ft[1]) + 12 * _sigma * (Fi[0]*Fi[0] + Fi[1]*Fi[1]) / (psi * psi * ((1/_Tau[1])-0.5));
	S[2] = -6 * (u[0]*Ft[0] + u[1]*Ft[1]) - 12 * _sigma * (Fi[0]*Fi[0] + Fi[1]*Fi[1]) / (psi * psi * ((1/_Tau[2])-0.5));
	S[3] = Ft[0];
	S[4] = -Ft[0];
	S[5] = Ft[1];
	S[6] = -Ft[1];
	S[7] = 2 * (u[0]*Ft[0] - u[1]*Ft[1]);
	S[8] = u[0]*Ft[1] + u[1]*Ft[0];




	// Collision in moment space
	
	for( uint k = 0 ; k < q ; k++ ) {

	    m[k] = m[k]  -  _Tau[k]*( m[k] - m_eq[k] )  +  ( 1 - 0.5*_Tau[k] ) * S[k];
	    
	}

	
	
	// Back to population space
	
	invM.matDotVec(m, _pdf[id]);		
	

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