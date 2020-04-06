#include <intForce.H>

using namespace std;


/** Force creator */

interactionForce* intForce::create (const string& dictName, const string& eqName, const latticeMesh& mesh, timeOptions& Time) {


    // Initialize types

    _ifMapType["singleRange"]                  = ifType::singleRange;
    _ifMapType["singleRangeMixed"]             = ifType::singleRangeMixed;
    _ifMapType["singleRangeWithContact"]       = ifType::singleRangeContact;
    _ifMapType["singleRangeMixedWithContact"]  = ifType::singleRangeMixedContact;

    
    
    // Load model name from dictionary

    dictionary dict(dictName);

    string itype = dict.lookUp<string>( eqName + "/Forces/Interaction/type" );




    // Assign interaction force
    
    if( _ifMapType.find(itype) != _ifMapType.end() ) {
	

	switch( _ifMapType[itype] ) {
	

	case ifType::singleRange:

	    return new singleRangeIntForce(dictName, eqName, mesh, Time);

	    break;


	case ifType::singleRangeMixed:

	    return new singleRangeMixedIntForce(dictName, eqName, mesh, Time);

	    break;


	case ifType::singleRangeContact:

	    return new singleRangeWithContact(dictName, eqName, mesh, Time);

	    break;


	case ifType::singleRangeMixedContact:

	    return new singleRangeMixedWithContact(dictName, eqName, mesh, Time);

	    break;	    	    
	    
	    
	}

    
    }

    
    else {
    
	cout << endl << " [ERROR]  Interaction type " << itype << " not available" << endl << endl;

	exit(1);

    }    



   
    
    return 0;


}
