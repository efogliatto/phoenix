#include <intForce.H>

using namespace std;


/** Force creator */

interactionForce* intForce::create (const string& dictName, const string& eqName, const latticeMesh& mesh, timeOptions& Time) {

    
    // Load model name from dictionary

    dictionary dict(dictName);

    string itype = dict.lookUp<string>( eqName + "/Forces/Interaction/type" );

    
    if( itype == "singleRange" ) {

    	return new singleRangeIntForce(dictName, eqName, mesh, Time);

    }

    else {

	if( itype == "singleRangeMixed" ) {

	    return new singleRangeMixedIntForce(dictName, eqName, mesh, Time);

	}


	else {


	    // Default
    
	    cout << endl << " [ERROR]  Interaction type " << itype << " not available" << endl << endl;

	    exit(1);


	}
    
    }

    
    return 0;


}
