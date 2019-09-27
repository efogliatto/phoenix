#include <adsForceCreator.H>

#include <dictionary.H>

using namespace std;


/** Force creator */
    
adhesiveForce* adsForceCreator::create( const string& dictName,
					const string& eqName,
					const latticeMesh& mesh,
					timeOptions& Time,
					const interactionForce* Fi) {


    // Initialize types

    _adsMapType["none"]           = adsType::NONE;
    _adsMapType["Phi-based-mod"]  = adsType::PHI_BASED_MOD;
    _adsMapType["randomSpots"]    = adsType::RND_SPOTS;    

    
    
    // Load model name from dictionary

    dictionary dict(dictName);

    string ftype = dict.lookUpOrDefault<string>( eqName + "/Forces/Adhesive/type", "none" );


    if( _adsMapType.find(ftype) != _adsMapType.end() ) {

	
    	switch( _adsMapType[ftype] ) {

	    
    	case adsType::NONE:

    	    return new noAds(dictName, eqName, mesh, Fi, Time);

    	    break;



    	case adsType::PHI_BASED_MOD:

    	    return new phiBasedMod(dictName, eqName, mesh, Fi, Time);

    	    break;	   


	    
    	case adsType::RND_SPOTS:

    	    return new rndSpots(dictName, eqName, mesh, Fi, Time);

    	    break;	   
	    
	    
    	}

	
    }

    
    else {


    	// Default
    
    	cout << endl << " [ERROR]  Adhesion force type " << ftype << " not available" << endl << endl;

    	exit(1);
    
    }

    

    
    return 0;


}
