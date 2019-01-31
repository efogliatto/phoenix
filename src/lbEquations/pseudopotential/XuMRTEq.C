#include <XuMRTEq.H>

using namespace std;


/** Default constructor */

XuMRTEq::XuMRTEq( const string& name,
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
    
    
    
    // Set equilibrium pdf values only if time is start time

    if( Time.currentTime() == Time.timeList()[0] )  {
	
	XuMRTEq::setEquilibrium();

    }
	

}
    

/** Default destructor */

XuMRTEq::~XuMRTEq() {}





/** Collision process */

const void XuMRTEq::collision() {

    
    // Lattice constants

    scalarMatrix M = mesh.lmodel()->MRTMatrix();

    scalarMatrix invM = mesh.lmodel()->MRTInvMatrix();    

    const uint nodes = mesh.npoints();

    const scalar cs2 = mesh.lmodel()->cs2();

    const scalar q = mesh.lmodel()->q();

    scalar Ft[3];

    
    // Partial distributions
    
    vector<scalar> m(15);     // m:  momentum space
    
    vector<scalar> m_eq(15);  // meq: equilibrium in momentum space

    vector<scalar> S(15);     // MRT force

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
	
	m_eq[0]  = r;     
	m_eq[1]  = r * (-1.0 + umag);
	m_eq[2]  = r * (1.0 - 5.0*umag);
	m_eq[3]  = r * u[0];
	m_eq[4]  = r * -7.0 * u[0] / 3.0;
	m_eq[5]  = r * u[1];
	m_eq[6]  = r * -7.0 * u[1] / 3.0;
	m_eq[7]  = r * u[2];
	m_eq[8]  = r * -7.0 * u[2] / 3.0;
	m_eq[9]  = r * ( 2.0 * u[0] * u[0]  -  u[1] * u[1]  -  u[2] * u[2] );
	m_eq[10] = r * ( u[1] * u[1]  -  u[2] * u[2] );
	m_eq[11] = r * u[0]*u[1];
	m_eq[12] = r * u[1]*u[2];
	m_eq[13] = r * u[0]*u[2];
	m_eq[14] = 0;

	

	// Distribution in moment space

	M.matDotVec( _pdf[id], m );


	
	// Source term

	scalar psi = F.potential( r, T.at(id), cs2 );

	F.total(Ft, id);

	const scalar Fi[3] = { F.interaction(id,0), F.interaction(id,1), F.interaction(id,2) };

	scalar uDotF = u[0]*Ft[0] + u[1]*Ft[1] + u[2]*Ft[2];

	scalar Fmag(0);

	for( uint j = 0 ; j < 3 ; j++ )	
	    Fmag += Fi[j] * Fi[j];

	

	S[0] = 0;
	S[1] = 2.0 * uDotF   +   (6.0 * _sigma * Fmag) / (psi*psi*((1/_Tau[1])-0.5));
	S[2] = -10.0 * uDotF;
	S[3] = Ft[0];
	S[4] = -7.0 * Ft[0] / 3.0;
	S[5] = Ft[1];
	S[6] = -7.0 * Ft[1] / 3.0;
	S[7] = Ft[2];
	S[8] = -7.0 * Ft[2] / 3.0;
	S[9] = 4.0 * u[0]*Ft[0]  -  2.0 * u[1]*Ft[1]  -  2.0 * u[2]*Ft[2];
	S[10] = 2.0 * u[1]*Ft[1]  -  2.0 * u[2]*Ft[2];
	S[11] = u[0]*Ft[1]  +  u[1]*Ft[0];
	S[12] = u[1]*Ft[2]  +  u[2]*Ft[1];
	S[13] = u[0]*Ft[2]  +  u[2]*Ft[0];
	S[14] = 0;
      




	// Collision in moment space
	
	for( uint k = 0 ; k < q ; k++ ) {

	    m[k] = m[k]  -  _Tau[k]*( m[k] - m_eq[k] )  +  ( 1 - 0.5*_Tau[k] ) * S[k];
	    
	}

	
	
	// Back to population space
	
	invM.matDotVec(m, _pdf[id]);		
	

    }

    

}




/** Set pdf to equilibrium values */

const void XuMRTEq::setEquilibrium() {

    
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
