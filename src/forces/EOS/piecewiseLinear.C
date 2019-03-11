#include <piecewiseLinear.H>

#include <dictionary.H>


using namespace std;


/** Default constructor */

piecewiseLinear::piecewiseLinear( const std::string& dictName,
				  const std::string& eqName )

    : EOS("piecewiseLinear") {


    // Open dictionary

    dictionary dict( dictName );

    _rho_1 = dict.lookUp<scalar>( eqName + "/Forces/EOS/rho_1" );

    _rho_2 = dict.lookUp<scalar>( eqName + "/Forces/EOS/rho_2" );

}


/** Default destructor */

piecewiseLinear::~piecewiseLinear() {}



/** EOS pressure */

const scalar piecewiseLinear::p_eos(const scalar& rho, const scalar& T) const {

    scalar p(0);

    if( rho < _rho_1 ) {

	p = rho * 0.64 / 3;
	
    }

    else {

	if(  ( rho > _rho_1 )  &&   ( rho < _rho_2 )){

	    p = _rho_1 * 0.64 / 3.0  -  (rho - _rho_1) * 0.04 / 3.0;
	
	}

	else {

	    p = _rho_1 * 0.64 / 3.0  -  (_rho_2 - _rho_1) * 0.04 / 3.0  +  (rho - _rho_2) / 3.0;

	}

    }
    

    return p;

}




/** Pressure derivative */

const scalar piecewiseLinear::dp_dT(const scalar& rho, const scalar& T) const {

    return 0;

}
