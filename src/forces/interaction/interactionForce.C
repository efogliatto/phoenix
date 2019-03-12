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


    // Read flag

    string onbnd = dict.lookUpOrDefault<string>( eqName + "/Forces/Interaction/OnBoundaries", "false" );

    if( onbnd == "true" )
    	_computeOnBnd = true;


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
