#include <fixedDensityBForce.H>


using namespace std;


/** Constructor */

fixedDensityBForce::fixedDensityBForce( const string& dictName,
					const string& eqName )
    : buoyantForce(dictName, eqName) {


    // Update reference density
    
    dictionary dict(dictName);

    _rhoRef = dict.lookUp<scalar>( eqName + "/Buoyancy/rhoRef" );    

}


/** Destructor */

fixedDensityBForce::~fixedDensityBForce() {}


/** Force */

const vector<scalar> fixedDensityBForce::force(const scalar& rho) {

    return { (rho-_rhoRef)*_gravity[0],
	     (rho-_rhoRef)*_gravity[1],
	     (rho-_rhoRef)*_gravity[2] };

}
