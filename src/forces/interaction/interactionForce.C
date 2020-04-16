#include <interactionForce.H>

#include <numeric>

using namespace std;



/** Constructor */

interactionForce::interactionForce( const string& dictName,
				    const string& eqName,
				    const latticeMesh& mesh,
				    timeOptions& Time )

    : _mesh(mesh),
      _force(mesh, Time, "Fi", IO::NO_READ, IO::NO_WRITE)// ,
      // _computeOnBnd(false)
{


    // Create eos

    EOSCreator creator;

    eos = creator.create(dictName, eqName);


    // Read Main interaction strength

    dictionary dict(dictName);

    _G = dict.lookUp<scalar>( eqName + "/Forces/Interaction/G" );


    // // Compute on boundaries

    // string onbnd = dict.lookUpOrDefault<string>( eqName + "/Forces/Interaction/OnBoundaries", "false" );

    // if( onbnd == "true" )
    // 	_computeOnBnd = true;


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


		hmin[bd.first] = dict.lookUpOrDefault<scalar>( eqname + "/Forces/Interaction/ContactAngle/" + bd.first + "/Hysteresis/min", 0) * M_PI / 180.0;

		hmax[bd.first] = dict.lookUpOrDefault<scalar>( eqname + "/Forces/Interaction/ContactAngle/" + bd.first + "/Hysteresis/max", 0) * M_PI / 180.0;		
		

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



/** Update contact nodes */

const scalar interactionForce::apparentContactAngle( const scalarField& rho, const std::string& bname ) {

    
    // Temporary contact nodes
    
    vector<scalar> cn;

    scalar _rhoAvgInt(3.484);

    const vector< vector<int> >& nb = _mesh.nbArray();
    
    
    
    // Move over all boundary nodes and detect density change

    const map< string, vector<uint> >& boundary = _mesh.boundaries();

    if( boundary.at(bname).size() > 0 ) {

	for( uint j = 0 ; j < boundary.at(bname).size()-1 ; j++ ) {

	    scalar y0( rho.at( boundary.at(bname)[j] ) - _rhoAvgInt ),
		   y1( rho.at( boundary.at(bname)[j+1] ) - _rhoAvgInt );

	    if( (y1 * y0)   <= 0) {

		int id = boundary.at(bname)[j];

		scalar angle = M_PI/2 - atan( -(rho.at(nb[id][4]) - rho.at(id)) / abs(rho.at(nb[id][3]) - rho.at(nb[id][1])) );
		
		cn.push_back( angle );


		
		id = boundary.at(bname)[j+1];

		angle = M_PI/2 - atan( - 2 * (rho.at(nb[id][4]) - rho.at(id)) / abs(rho.at(nb[id][3]) - rho.at(nb[id][1])) );
		
		cn.push_back( angle );
		
	    }

	}

    }

    

    return accumulate(cn.begin(), cn.end(), 0.0) / cn.size();
 

}
