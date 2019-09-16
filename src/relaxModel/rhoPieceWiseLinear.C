#include <rhoPieceWiseLinear.H>

#include <dictionary.H>

using namespace std;


/** Default constructor */

rhoPieceWiseLinear::rhoPieceWiseLinear( const std::string& entry ) : relaxModel(entry) {

    _name = "uniform";

    
    dictionary dict("properties/macroProperties");

    _Tau_v = dict.lookUp< vector<scalar> >( entry + "/Tau_v" );

    _Tau_l = dict.lookUp< vector<scalar> >( entry + "/Tau_l" );
    

    _rho_v = dict.lookUp<scalar>( entry + "/rho_v" );

    _rho_l = dict.lookUp<scalar>( entry + "/rho_l" );
    
}


/** Destructor */

rhoPieceWiseLinear::~rhoPieceWiseLinear() {}



/** Density-dependent relaxation factor */

const scalar rhoPieceWiseLinear::tau( const scalar rho, const uint i = 0 ) {

    scalar t(0);
    
    if (rho <= _rho_v) {

	t = _Tau_v[i];
	
    }

    else {

	if (rho >= _rho_l) {

	    t = _Tau_l[i];
	
	}

	else {

	    t = (_Tau_l[i] - _Tau_v[i]) * (rho - _rho_v) / (_rho_l - _rho_v) + _Tau_v[i];

	}

    }



    return t;
	
}
