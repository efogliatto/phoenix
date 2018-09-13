#include <EOSCreator.H>

#include <iostream>

#include <dictionary.H>


using namespace std;


EOS* EOSCreator::create( const string& dictName, const string& eqName ) {

    
    // Load model name from dictionary

    dictionary dict(dictName);

    string etype = dict.lookUp<string>( eqName + "/Forces/EOS/type" );

    
    if( etype == "vanDerWaals" ) {

    	return new vanDerWaals(dictName, eqName);

    }

    else {


    	// Default
    
    	cout << endl << " [ERROR]  EOS " << etype << " not available" << endl << endl;

    	exit(1);
    
    }

    
    return 0;

}
