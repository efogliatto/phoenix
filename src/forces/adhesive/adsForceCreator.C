#include <adsForceCreator.H>

#include <dictionary.H>

using namespace std;


/** Force creator */
    
adhesiveForce* adsForceCreator::create( const string& dictName,
					const string& eqName,
					const latticeMesh& mesh ) {

    
    // Load model name from dictionary

    dictionary dict(dictName);

    string ftype = dict.lookUpOrDefault<string>( eqName + "/Forces/Adhesive/type", "none" );

    
    if( ftype == "liAdhesive" ) {

    	return new liAdhesive(dictName, eqName, mesh);

    }

    else {

	if( ftype == "none" ) {

	    return new noAds(dictName, eqName, mesh);

	}


	else {


	    // Default
    
	    cout << endl << " [ERROR]  Adhesive force type " << ftype << " not available" << endl << endl;

	    exit(1);


	}
    
    }

    
    return 0;


}
