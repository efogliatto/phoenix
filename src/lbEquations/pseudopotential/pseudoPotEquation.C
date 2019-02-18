#include <pseudoPotEquation.H>

using namespace std;


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
      F("properties/macroProperties", name, mesh_, Time_, rho_, T_) {


    
    // Update Forces

    F.update(rho,T);



    // // Read Boundary conditions

    // dictionary dict("start/boundaries");

    // const map< string, vector<uint> >& bnd = mesh.boundaries();

    // ppBndCreator BndCreator;

    // for(map< string, vector<uint> >::const_iterator iter = bnd.begin(); iter != bnd.end(); ++iter)  {
	
    // 	string bdname = iter->first;

    // 	_boundaries.push_back( BndCreator.create(name, bdname, mesh.boundaryNodes(bdname) ) );

    // }

}


/** Default destructor */

pseudoPotEquation::~pseudoPotEquation() {}










/** Compute local density */

const scalar pseudoPotEquation::localDensity( const uint& id ) const {

    scalar r(0);

    const uint q = mesh.lmodel()->q();
    
    for( uint k = 0 ; k < q ; k++ )
    	r += _pdf[id][k];


    return r;    

}


/** Compute local velocity */

const void pseudoPotEquation::localVelocity( vector<scalar>& v, const uint& id, const bool updDens ) const {
    

    // Local velocity
    
    vector<scalar> lv = {0,0,0};


    // Local density

    scalar localRho(0);

    if(updDens) {

	localRho = localDensity(id);

    }

    else {

	localRho = rho.at(id);

    }


    // Lattice constants

    const vector< vector<int> >& vel = mesh.lmodel()->lvel();

    const uint q = mesh.lmodel()->q();


    // Compute first moment
    
    for( uint j = 0 ; j < 3 ; j++ ) {

    	for( uint k = 0 ; k < q ; k++ ) {

    	    lv[j] += vel[k][j] * _pdf[id][k];
		    
    	}
	    
    }


    // Add interaction force and divide by density

    scalar Ft[3];

    F.total(Ft, id);
    
    for( uint j = 0 ; j < 3 ; j++ )
    	lv[j] = ( lv[j]   +   0.5 * Ft[j]   ) / localRho;
	


    // Copy to global array
    for( uint j = 0 ; j < 3 ; j++ ) 
    	v[j] = lv[j];
	

}







/** Collision process */

const void pseudoPotEquation::collision() {}



/** Update macroscopic density */

const void pseudoPotEquation::updateMacroDensity() {

    for( uint i = 0 ; i < mesh.npoints() ; i++ )
	rho[i] = localDensity(i);

}


const void pseudoPotEquation::updateMacroDensity(const uint& first, const uint& last) {

    for( uint i = first ; i < last ; i++ )
	rho[i] = localDensity(i);

}




/** Update macroscopic velocity */
/** Repeat for efficiency */

const void pseudoPotEquation::updateMacroVelocity() {

    
    // // Update forces first

    // F.update(rho,T);



    // Lattice constants

    const vector< vector<int> > vel = mesh.lmodel()->lvel();

    const scalar q = mesh.lmodel()->q();

    
    // Local velocity
    
    scalar lv[3] = {0,0,0};

    scalar Ft[3];
    

    
    for( uint i = 0 ; i < mesh.npoints() ; i++ ) {


	// Compute first moment
    
	for( uint j = 0 ; j < 3 ; j++ ) {

	    lv[j] = 0;
	    
	    for( uint k = 0 ; k < q ; k++ ) {

		lv[j] += vel[k][j] * _pdf[i][k];
		    
	    }
	    
	}


	// Add interaction force and divide by density

	F.total(Ft, i);
    
	for( uint j = 0 ; j < 3 ; j++ )
	    U[i][j] = ( lv[j]   +   0.5 * Ft[j]   ) / rho.at(i);

	
    }

}





/** Forced interface for equilibrium distribution */

const void pseudoPotEquation::eqPS( std::vector<scalar>& n, const scalar& rho_, const std::vector<scalar>& U_ ) const {

    
    // Lattice constants
    
    const uint q = mesh.lmodel()->q();

    const scalar cs2 = mesh.lmodel()->cs2();

    const vector<scalar>& omega = mesh.lmodel()->omega();

    const vector< vector<int> >& vel = mesh.lmodel()->lvel();


    scalar umag(0);

    for( uint j = 0 ; j < 3 ; j++ )
	umag += U_[j] * U_[j];
    
    
    for( uint k = 0 ; k < q ; k++ ) {

	// Dot product

	scalar alpha(0);
	
	for( uint j = 0 ; j < 3 ; j++ )
	    alpha += vel[k][j] * U_[j];

	    
	n[k] = rho_ * omega[k] * ( 1 + alpha/cs2   +   0.5 * alpha * alpha / (cs2*cs2)  -  0.5 * umag / cs2 );

    }

}



/** Update potential as scalar field */

const void pseudoPotEquation::updatePotential( scalarField& phi ) {

    const scalar cs2 = mesh.lmodel()->cs2();
    
    for( uint i = 0 ; i < mesh.npoints() ; i++ ) {

	phi[i] = F.potential( rho.at(i), T.at(i), cs2 );

    }

}



/** Compute and set pressure field */

const void pseudoPotEquation::pressure( const scalarField& phi, scalarField& p ) {

    // Base case. Ideal LBM

    const scalar cs2 = mesh.lmodel()->cs2();
    
    for( uint i = 0 ; i < mesh.npoints() ; i++ ) {

	p[i] = rho.at(i) * cs2;

    }    

}
