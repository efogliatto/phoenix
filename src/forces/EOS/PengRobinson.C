#include <PengRobinson.H>

#include <dictionary.H>


using namespace std;


/** Default constructor */

PengRobinson::PengRobinson( const std::string& dictName,
			    const std::string& eqName )

    : EOS("PengRobinson") {


    // Open dictionary

    dictionary dict( dictName );

    a = dict.lookUp<scalar>( eqName + "/Forces/EOS/a" );

    b = dict.lookUp<scalar>( eqName + "/Forces/EOS/b" );

    w = dict.lookUp<scalar>( eqName + "/Forces/EOS/w" );    

}


/** Default destructor */

PengRobinson::~PengRobinson() {}


/** EOS pressure */

const scalar PengRobinson::p_eos(const scalar& rho, const scalar& T) const {


    scalar Tc = 0.0778 * a / (0.45724 * b);

    scalar theta = 1 + (0.37464 + 1.54226 * w - 0.26992 * w * w) * (1 - sqrt(T/Tc));

    theta = theta * theta;

    

    scalar p = rho * T / (1.0 - rho * b)  -  a * theta * rho * rho / ( 1  +  2 * b * rho  -  b * b * rho * rho );
    
    return  p;

}



/** Pressure derivative */

const scalar PengRobinson::dp_dT(const scalar& rho, const scalar& T) const {

    scalar Tc = 0.0778 * a / (0.45724 * b);

    scalar theta = 1 + (0.37464 + 1.54226 * w - 0.26992 * w * w);

    theta = theta * theta * (1-sqrt(Tc/T)) / Tc;
    

    scalar p = rho / (1.0 - rho * b) - a * theta * rho * rho / ( 1  +  2 * b * rho  -  b * b * rho * rho );
    
    return p;

}
