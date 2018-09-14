#include <externalForce.H>

using namespace std;


/** Constructor */

externalForce::externalForce( const string& dictName, const string& eqName ) {

    // Update force value
    
    dictionary dict(dictName);

    _force = dict.lookUp< vector<scalar> >( eqName + "/External/value" );


}


/** Destructor */

externalForce::~externalForce() {}
