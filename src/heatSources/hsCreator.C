#include <hsCreator.H>

using namespace std;


/** Force creator */

heatSource* hsCreator::create (const string& dictName, const string& eqName, const latticeMesh& mesh, timeOptions& Time) {

    
    // Load model name from dictionary

    dictionary dict(dictName);

    string itype = dict.lookUp<string>( eqName + "/HeatSource/type" );

    
    if( itype == "Markus-Hazi" ) {

    	return new markusHaziHS(dictName, eqName, mesh, Time);

    }

    else {

	// Default
    
	cout << endl << " [ERROR]  Heat Source type " << itype << " not available" << endl << endl;

	exit(1);

    
    }

    
    return 0;


}
