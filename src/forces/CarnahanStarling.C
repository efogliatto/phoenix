#include <CarnahanStarling.H>

#include <dictionary.H>


using namespace std;


/** Default constructor */

CarnahanStarling::CarnahanStarling( const std::string& dictName,
			  const std::string& eqName )

    : EOS("CarnahanStarling") {


    // Open dictionary

    dictionary dict( dictName );

    a = dict.lookUp<scalar>( eqName + "/Forces/EOS/a" );

    b = dict.lookUp<scalar>( eqName + "/Forces/EOS/b" );

}


/** Default destructor */

CarnahanStarling::~CarnahanStarling() {}


/** EOS pressure */

const scalar CarnahanStarling::p_eos(const scalar& rho, const scalar& T) const {


    scalar alpha = b * rho / 4.0;

    scalar beta  = 1 - alpha;

    scalar p = rho * T * ( (1 + alpha + alpha*alpha - alpha*alpha*alpha) / ( beta*beta*beta )  ) - a * rho * rho;

    
    return  p;

}



/** Pressure derivative */

const scalar CarnahanStarling::dp_dT(const scalar& rho, const scalar& T) const {

    return rho / (1 - rho * b);

}
