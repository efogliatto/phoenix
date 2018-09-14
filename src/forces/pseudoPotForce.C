#include <pseudoPotForce.H>

using namespace std;


/** Constructor */

pseudoPotForce::pseudoPotForce( const string& dictName, const string& eqName, const latticeMesh& mesh, timeOptions& Time ) {


    // Create interaction force

    intForce ifcreator;

    _Fi = ifcreator.create(dictName, eqName, mesh, Time);


    
    // Create buoyant force

    bForces bfcreator;

    _Fb = bfcreator.create(dictName, eqName);
    

}



/** Destructor */

pseudoPotForce::~pseudoPotForce() {}



/** Total force  at node*/

const vector<scalar> pseudoPotForce::total( const uint& id ) const {

    vector<scalar> Fi = _Fi->force(id);

    vector<scalar> Fb = _Fb->force(id);

    return { Fi[0] + Fb[0],
	     Fi[1] + Fb[1],
	     Fi[2] + Fb[2]};    

}



/** Interaction force at node */

const vector<scalar> pseudoPotForce::interaction( const uint& id ) const {

    return _Fi->force(id);

}



/** Interaction force at node */

const vector<scalar> pseudoPotForce::buoyant( const uint& id ) const {

    return _Fb->force(id);

}



/** Update forces */

void pseudoPotForce::update( scalarField& rho, scalarField& T ) {

    _Fi->update(rho,T);

}



/** Sync force fields */

void pseudoPotForce::sync() {

    _Fi->sync();

}
