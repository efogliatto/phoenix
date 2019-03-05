#include <stCreator.H>

#include <dictionary.H>

using namespace std;


/** Surface Tension creator */
    
surfaceTension* stCreator::create( const string& dictName,
				   const string& eqName,
				   const latticeMesh& mesh ) {

    
    // Load model name from dictionary

    dictionary dict(dictName);

    string sttype = dict.lookUpOrDefault<string>( eqName + "/Forces/SurfaceTension/type", "none" );

    
    if( sttype == "none" ) {

    	return new noSurfaceTension(dictName, eqName, mesh);

    }

    else {

	// if( ftype == "none" ) {

	//     return new noAds(dictName, eqName, mesh);

	// }

	// else {


	//     if( ftype == "simpleAdhesive" ) {

	// 	return new simpleAdhesive(dictName, eqName, mesh);

	//     }

	    
	//     else {


		// Default
    
		cout << endl << " [ERROR]  Surface tension type " << sttype << " not available" << endl << endl;

		exit(1);


	//     }

	// }
    
    }

    
    return 0;


}
