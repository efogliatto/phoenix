#include <uniformTau.H>

#include <dictionary.H>

using namespace std;


/** Default constructor */

uniformTau::uniformTau( const std::string& entry ) : relaxModel(entry) {

    _name = "uniform";

    
    dictionary dict("properties/macroProperties");

    _Tau = dict.lookUp< vector<scalar> >( entry + "/Tau" );

}


/** Destructor */

uniformTau::~uniformTau() {}


/** Density-dependent relaxation factor */

const scalar uniformTau::tau( const scalar rho, const uint i = 0 ) {

    return _Tau[i];

}
