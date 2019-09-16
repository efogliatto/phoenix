#include <rhoPieceWise.H>

#include <dictionary.H>

using namespace std;


/** Default constructor */

rhoPieceWise::rhoPieceWise( const std::string& entry ) : relaxModel(entry) {

    _name = "rhoPieceWise";

    
    dictionary dict("properties/macroProperties");

    _Tau_v = dict.lookUp< vector<scalar> >( entry + "/Tau_v" );

    _Tau_l = dict.lookUp< vector<scalar> >( entry + "/Tau_l" );
    

    _rho_lim = dict.lookUp<scalar>( entry + "/rho_lim" );
    
}


/** Destructor */

rhoPieceWise::~rhoPieceWise() {}



/** Density-dependent relaxation factor */

const scalar rhoPieceWise::tau( const scalar rho, const uint i = 0 ) {

    scalar t(0);
    
    if (rho <= _rho_lim) {

	t = _Tau_v[i];
	
    }

    else {

	t = _Tau_l[i];	

    }


    return t;
	
}
