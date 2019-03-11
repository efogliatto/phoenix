#include <pseudoPotForce.H>

using namespace std;


/** Constructor */

pseudoPotForce::pseudoPotForce( const string& dictName,
				const string& eqName,
				const latticeMesh& mesh,
				timeOptions& Time,
				const scalarField& rho,
				const scalarField& T)

    : _Fe(dictName, eqName),
      _rho(rho),
      _T(T),
      _mesh(mesh){


    // Create interaction force

    intForce ifcreator;

    _Fi = ifcreator.create(dictName, eqName, mesh, Time);


    
    // Create buoyant force

    bForces bfcreator;

    _Fb = bfcreator.create(dictName, eqName);


    
    // Create adhesive force

    adsForceCreator adcreator;

    _Fads = adcreator.create(dictName, eqName, mesh);



    // Create additional surface tension term

    stCreator st;

    _St = st.create(dictName, eqName, mesh);
    

}



/** Destructor */

pseudoPotForce::~pseudoPotForce() {}



/** Total force  at node*/

const vector<scalar> pseudoPotForce::total( const uint& id ) const {


    // // Interaction
    
    // const scalar Fi[3] = { _Fi->force(id,0), _Fi->force(id,1), _Fi->force(id,2) };


    // // Bouyant
    
    // const scalar Fb[3] = { _Fb->force( _rho.at(id), 0 ), _Fb->force( _rho.at(id), 1 ), _Fb->force( _rho.at(id), 2 ) };


    // // Adhesion

    // scalar Fads[3] = {0,0,0};

    // _Fads->force(id, _rho.at(id), _T.at(id), Fads);

    
    
    // return { Fi[0] + Fb[0] + Fads[0] + _Fe[0] * _rho.at(id),
    // 	     Fi[1] + Fb[1] + Fads[1] + _Fe[1] * _rho.at(id),
    // 	     Fi[2] + Fb[2] + Fads[2] + _Fe[2] * _rho.at(id) };


    scalar F[3] = {0,0,0};

    pseudoPotForce::total(F,id);

    return {F[0], F[1], F[2]};
    

}



/** Total force  at node*/

void pseudoPotForce::total( scalar Ft[3], const uint& id ) const {

    // if( _mesh.latticePoint(id)[1] != 0 ) {


	// Interaction
    
	const scalar Fi[3] = { _Fi->force(id,0), _Fi->force(id,1), _Fi->force(id,2) };


	// Buoyant
    
	const scalar Fb[3] = { _Fb->force( _rho.at(id), 0 ), _Fb->force( _rho.at(id), 1 ), _Fb->force( _rho.at(id), 2 ) }; 

    
	// Adhesion

	scalar Fads[3] = {0,0,0};

	_Fads->force(id, _rho, _T, Fads);

    
    
	Ft[0] = Fi[0] + Fb[0] + Fads[0] + _Fe[0] * _rho.at(id);
	Ft[1] = Fi[1] + Fb[1] + Fads[1] + _Fe[1] * _rho.at(id);
	Ft[2] = Fi[2] + Fb[2] + Fads[2] + _Fe[2] * _rho.at(id);

    // }

    // else {

    // 	Ft[0] = 0;
    // 	Ft[1] = 0;
    // 	Ft[2] = 0;	

    // }

   
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

    _Fb->update(rho);

}



/** Sync force fields */

void pseudoPotForce::sync() {

    _Fi->sync();

}



/** Interaction potential */

const scalar pseudoPotForce::potential(const scalar& rho, const scalar& T, const scalar& cs2) const {

    return _Fi->potential(rho,T,cs2);

}



/** Adequately signed potential strength */

const scalar pseudoPotForce::signedPotential( const scalar& rho, const scalar& T, const scalar& cs2 ) const {

    return _Fi->signedPotentialStrength(rho,T,cs2);

}



/** Additional surface tension term at node */

const void pseudoPotForce::addSurfaceTension( const uint& i, vector<scalar>& C, const std::vector<scalar>& Tau ) const {

    _St->ST( i, _rho, _T, C, _Fi, Tau );

}
