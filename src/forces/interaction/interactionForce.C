#include <interactionForce.H>

using namespace std;



/** Constructor */

interactionForce::interactionForce( const string& dictName,
				    const string& eqName,
				    const latticeMesh& mesh,
				    timeOptions& Time )

    : _mesh(mesh),
      _force(mesh, Time, "Fi", IO::NO_READ, IO::NO_WRITE),
      _computeOnBnd(false) {


    // Create eos

    EOSCreator creator;

    eos = creator.create(dictName, eqName);


    // Read Main interaction strength

    dictionary dict(dictName);

    _G = dict.lookUp<scalar>( eqName + "/Forces/Interaction/G" );


    // Compute on boundaries

    string onbnd = dict.lookUpOrDefault<string>( eqName + "/Forces/Interaction/OnBoundaries", "false" );

    if( onbnd == "true" )
    	_computeOnBnd = true;


    // Contact angle options

    readGeometricContact(eqName);


}



/** Destructor */

interactionForce::~interactionForce() {}



/** Interaction potential */

const scalar interactionForce::potential( const scalar& rho, const scalar& T, const scalar& cs2 ) const {

    scalar a = 2 * (eos->p_eos(rho,T) - rho * cs2)  /  _G;

    scalar b(0);

    (a >= 0)  ?	 b = sqrt(a)  :  b = sqrt(-a);
    
    return b;

}



/** Adecuately signed potential strength */

const scalar interactionForce::signedPotentialStrength( const scalar& rho, const scalar& T, const scalar& cs2 ) const {

    scalar a = 2 * (eos->p_eos(rho,T) - rho * cs2)  /  _G;
    
    scalar b(0);

    (a >= 0)  ?	 b = _G  :  b = -_G;
    
    return b;
    
}



/** Set force at specific node */

const void interactionForce::set( const uint& i, const vector<scalar>& Fint ) {

    for( uint j = 0 ; j < 3 ; j++ )
	_force[i][j] = Fint[j];

}



/** Read geometric contact angle properties */

const void interactionForce::readGeometricContact( const string& eqname ) {

    
    // Read contact angles for each boundary

    dictionary dict("properties/macroProperties");

    string contact = dict.lookUpOrDefault<string>( eqname + "/Forces/Interaction/ContactAngle/type", "none");
    

    if( contact == "geometric" ) {


	// Turn on flag

	_withGeomContact = true;

	

	// Read angles for each boundary

	const map< string, vector<uint> >& boundary = _mesh.boundaries();

	map<string, scalar> angle;

	map<string, scalar> hmin;

	map<string, scalar> hmax;
	

	for( const auto& bd : boundary) {

	    scalar ang = dict.lookUpOrDefault<scalar>( eqname + "/Forces/Interaction/ContactAngle/" + bd.first + "/Theta/", -1);

	    if(ang >= 0) {
		
		angle[bd.first] = ang * M_PI / 180.0;


		hmin[bd.first] = dict.lookUpOrDefault<scalar>( eqname + "/Forces/Interaction/ContactAngle/" + bd.first + "/Hysteresis/min", 0);

		hmax[bd.first] = dict.lookUpOrDefault<scalar>( eqname + "/Forces/Interaction/ContactAngle/" + bd.first + "/Hysteresis/max", 0);		
		

	    }

	}



	// Create map for each boundary node

	for( const auto& bd : angle) {

	    for( const auto& id : boundary.at(bd.first) ) {

		_contactAngle[id] = bd.second;

		_hysteresis[id] = {hmin[bd.first], hmax[bd.first]};

	    }

	}
	

    }

}
