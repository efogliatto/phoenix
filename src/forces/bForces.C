#include <bForces.H>

using namespace std;


buoyantForce* bForces::create( const string& dictName, const string& eqName ) {

    
    // Load model name from dictionary

    dictionary dict(dictName);

    string ftype = dict.lookUp<string>( eqName + "/Forces/Buoyancy/type" );

    
    if( ftype == "fixedDensity" ) {

    	return new fixedDensityBForce(dictName, eqName);

    }

    else {


    	// Default
    
    	cout << endl << " [ERROR]  Interaction type " << ftype << " not available" << endl << endl;

    	exit(1);
    
    }

    
    return 0;


}
