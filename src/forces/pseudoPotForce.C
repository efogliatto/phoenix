#include <pseudoPotForce.H>

using namespace std;


/** Constructor */

pseudoPotForce::pseudoPotForce( const string& dictName,
				const string& eqName,
				const latticeMesh& mesh,
				timeOptions& Time,
				const scalarField& rho)

    : _Fe(dictName, eqName),
      _rho(rho) {


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

    const scalar Fi[3] = { _Fi->force(id,0), _Fi->force(id,1), _Fi->force(id,2) };

    const scalar Fb[3] = { _Fb->force( _rho.at(id), 0 ), _Fb->force( _rho.at(id), 1 ), _Fb->force( _rho.at(id), 2 ) }; 
    
    return { Fi[0] + Fb[0] + _Fe[0],
	     Fi[1] + Fb[1] + _Fe[1],
	     Fi[2] + Fb[2] + _Fe[2] };

}



/** Total force  at node*/

void pseudoPotForce::total( scalar Ft[3], const uint& id ) {

    const scalar Fi[3] = { _Fi->force(id,0), _Fi->force(id,1), _Fi->force(id,2) };

    const scalar Fb[3] = { _Fb->force( _rho.at(id), 0 ), _Fb->force( _rho.at(id), 1 ), _Fb->force( _rho.at(id), 2 ) }; 
    
    Ft[0] = Fi[0] + Fb[0] + _Fe[0];
    Ft[1] = Fi[1] + Fb[1] + _Fe[1];
    Ft[2] = Fi[2] + Fb[2] + _Fe[2];

}



/** Interaction force at node */

const vector<scalar>& pseudoPotForce::interaction( const uint& id ) const {

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



/** Interaction potential */

const scalar pseudoPotForce::potential(const scalar& rho, const scalar& T, const scalar& cs2) const {

    return _Fi->potential(rho,T,cs2);

}
