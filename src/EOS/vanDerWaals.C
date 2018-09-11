#include <vanDerWaals.H>

#include <dictionary.H>


using namespace std;


/** Default constructor */

vanDerWaals::vanDerWaals( const std::string fname ) : EOS("vanDerWaals") {


    // Open dictionary

    dictionary dict( fname );

    a = dict.lookUp<scalar>("EOS/a_vdw");

    b = dict.lookUp<scalar>("EOS/b_vdw");

}


/** Default destructor */

vanDerWaals::~vanDerWaals() {}


/** EOS pressure */

const scalar vanDerWaals::p_eos(const scalar& rho, const scalar& T) const {

    return  rho * T / (1 - rho * b)  -  a * rho * rho;

}
