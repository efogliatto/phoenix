#include <vanDerWaals.H>

#include <dictionary.H>


using namespace std;


/** Default constructor */

vanDerWaals::vanDerWaals( const std::string& dictName,
			  const std::string& eqName )

    : EOS("vanDerWaals") {


    // Open dictionary

    dictionary dict( dictName );

    a = dict.lookUp<scalar>( eqName + "/Forces/EOS/a" );

    b = dict.lookUp<scalar>( eqName + "/Forces/EOS/b" );

}


/** Default destructor */

vanDerWaals::~vanDerWaals() {}


/** EOS pressure */

const scalar vanDerWaals::p_eos(const scalar& rho, const scalar& T) const {

    return  rho * T / (1 - rho * b)  -  a * rho * rho;

}



/** Pressure derivative */

const scalar vanDerWaals::dp_dT(const scalar& rho, const scalar& T) const {

    return rho / (1 - rho * b);

}
