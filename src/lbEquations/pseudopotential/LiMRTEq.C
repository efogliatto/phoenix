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
    
    
    
    // Set equilibrium pdf values only if time is start time

    if( Time.currentTime() == Time.timeList()[0] )  {

	LiMRTEq::setEquilibrium();

    }

}
    

/** Default destructor */

LiMRTEq::~LiMRTEq() {}




/** Collision process */

const void LiMRTEq::collision() {

    
    // Lattice constants

    const scalarMatrix& M = mesh.lmodel()->MRTMatrix();

    const scalarMatrix& invM = mesh.lmodel()->MRTInvMatrix();    

    const scalar cs2 = mesh.lmodel()->cs2();

    const scalar q = mesh.lmodel()->q();

    scalar Ft[3];



    // Collision range depends on surface tension model
    
    uint nodes( mesh.npoints() );

    if( F.sTModel() == surfaceTension::stType::liST )
	nodes = mesh.local();

    

    
    // Partial distributions
    
    vector<scalar> m(9);     // m:  momentum space
    
    vector<scalar> m_eq(9);  // meq: equilibrium in momentum space

    vector<scalar> S(9);     // MRT force

    vector<scalar> C(9);     // Surface tension term


    
    // Local copy of relaxation factors

    vector<scalar> localTau(q);
    


    

    // Move over all points

    for( uint id = 0 ; id < nodes ; id++ ) {


	// Local velocity and density

	const scalar r = rho.at(id);

	const scalar u[3] = { U.at(id,0), U.at(id,1), U.at(id,2) };


	
	// Update local values of relaxation factors

	for(uint k = 0 ; k < q ; k++)
	    localTau[k] = _relax->tau(r,k);

	

	// Velocity magnitude

	scalar umag (0);
	
	for( uint j = 0 ; j < 3 ; j++ )	
	    umag += u[j] * u[j];


	
	// Compute equilibrium in moment space
	
	m_eq[0] = r;
	m_eq[1] = r * (-2.0 + 3.0*umag);
	m_eq[2] = r * (1.0 - 3.0*umag);
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


	S[0] = 0.0;
	S[1] =  6.0 * (u[0]*Ft[0] + u[1]*Ft[1]) + 12.0 * _sigma * (Fi[0]*Fi[0] + Fi[1]*Fi[1]) / (psi * psi * ((1.0/localTau[1])-0.5));
	S[2] = -6.0 * (u[0]*Ft[0] + u[1]*Ft[1]) - 12.0 * _sigma * (Fi[0]*Fi[0] + Fi[1]*Fi[1]) / (psi * psi * ((1.0/localTau[2])-0.5));
	S[3] = Ft[0];
	S[4] = -Ft[0];
	S[5] = Ft[1];
	S[6] = -Ft[1];
	S[7] = 2.0 * (u[0]*Ft[0] - u[1]*Ft[1]);
	S[8] = u[0]*Ft[1] + u[1]*Ft[0];



	
	// Surface tension extra term

	F.addSurfaceTension(id, C, localTau);
	


	
	// Collision in moment space
	
	for( uint k = 0 ; k < q ; k++ ) {

	    m[k] = m[k]
		-  localTau[k]*( m[k] - m_eq[k] )
		+  ( 1.0 - 0.5*localTau[k] ) * S[k] 
		+  C[k]
		;
	    
	}

	

	
	// Back to population space
	
	invM.matDotVec(m, _pdf[id]);		
	

    }


    
    // Extra sync only for Li's surface tension model
    
    if( F.sTModel() == surfaceTension::stType::liST )
	_pdf.sync();
	

    

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



/** Compute and update isotropic pressure as scalar field */

const void LiMRTEq::pressure( const scalarField& phi, scalarField& p ) {

    
    
    // Lattice constants   
    
    const scalar cs2 = mesh.lmodel()->cs2();

    scalar G(-1.0);

    const scalar c(1);


    // Potential gradient
    
    scalar grad[3] = {0,0,0};

    scalar phiMag(0);

    


    // Compute only at local points

    for(uint i = 0 ; i < mesh.local() ; i++){


	// Compute potential gradient

	phi.grad(grad, i);

	phiMag = 0;

	for( uint j = 0 ; j < 3 ; j++ )
	    phiMag += grad[j] * grad[j];


	// Signed potential strength

	G = F.signedPotential(rho.at(i), T.at(i), cs2);
	
	
	// Isotropic pressure
	
	p[i] = rho.at(i)*cs2
	    + (G * c * c / 2.0) * phi.at(i) * phi.at(i)
	    + (G * c * c * c * c/ 12.0) * phi.at(i) * phi.laplacian(i)
	    + 2 * G * G * c * c * c * c * _sigma * phiMag
	    ;

    }


    // Sync pressure field
    
    p.sync();
    
    

}
