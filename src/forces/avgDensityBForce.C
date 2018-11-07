#include <avgDensityBForce.H>

using namespace std;


/** Constructor */

avgDensityBForce::avgDensityBForce( const string& dictName,
				    const string& eqName )
    
    : buoyantForce(dictName, eqName) {


    // Update reference density
    


}


/** Destructor */

avgDensityBForce::~avgDensityBForce() {}


/** Force */

const vector<scalar> avgDensityBForce::force(const scalar& rho) {

    // return { (rho-_rhoRef)*_gravity[0],
    // 	     (rho-_rhoRef)*_gravity[1],
    // 	     (rho-_rhoRef)*_gravity[2] };

    
    const scalar a( rho*(1-_rhoRef/rho) );
    
    return { a*_gravity[0],
	     a*_gravity[1],
	     a*_gravity[2] };    

}



/** Update coefficients */

const void avgDensityBForce::update(const scalarField& rho) {

    _rhoRef = rho.average();

}
