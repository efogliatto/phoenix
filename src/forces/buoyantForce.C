#include <buoyantForce.H>


using namespace std;



/** Constructor */

buoyantForce::buoyantForce( const string& dictName, const string& eqName ) {


    // Update gravity value
    
    dictionary dict(dictName);

    _gravity = dict.lookUp< vector<scalar> >( eqName + "/Buoyancy/gravity" );
    

}


/** Destructor */

buoyantForce::~buoyantForce() {}
