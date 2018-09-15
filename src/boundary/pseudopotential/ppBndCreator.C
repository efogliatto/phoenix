#include <ppBndCreator.H>

using namespace std;


ppBndCond* ppBndCreator::create( const string& eqName, const string& bdName, const vector<uint>& nodes ) {

    
    // Load model name from dictionary

    dictionary dict("start/boundaries");

    string btype = dict.lookUp<string>( eqName + "/" + bdName + "/type" );

    
    if( btype == "fixedU" ) {

    	return new ppFixedU( eqName, bdName, nodes );

    }

    

    else {


	if( btype == "periodic" ) {

	    return new ppPeriodic( nodes );

	}

	else {


	    // Default
    
	    cout << endl << " [ERROR]  Boundary condition " << btype << " not available for pseudopotential model" << endl << endl;

	    exit(1);

	}
    
    }

    
    return 0;


}
