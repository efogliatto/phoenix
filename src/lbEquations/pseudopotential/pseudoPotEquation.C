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
      F("properties/macroProperties", name, mesh_, Time_, rho_){


    
    // Update Forces

    F.update(rho,T);



    // Read Boundary conditions

    dictionary dict("start/boundaries");

    const map< string, vector<uint> >& bnd = mesh.boundaries();

    ppBndCreator BndCreator;

    for(map< string, vector<uint> >::const_iterator iter = bnd.begin(); iter != bnd.end(); ++iter)  {
	
	string bdname = iter->first;

	_boundaries.push_back( BndCreator.create(name, bdname, mesh.boundaryNodes(bdname) ) );

    }

}


/** Default destructor */

pseudoPotEquation::~pseudoPotEquation() {}










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

const void pseudoPotEquation::localVelocity( vector<scalar>& v, const uint& id ) {
    

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

    scalar Ft[3];

    F.total(Ft, id);
    
    for( uint j = 0 ; j < 3 ; j++ )
    	lv[j] = ( lv[j]   +   0.5 * Ft[j]   ) / rho.at(id);
	


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

    
    // Update forces first

    F.update(rho,T);



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




/** Update boundaries */

void pseudoPotEquation::updateBoundaries() {

    for(uint i = 0 ; i < _boundaries.size() ; i++) {

	_boundaries[i]->update( mesh, _pdf, rho, T, U, &F );

    }

}
