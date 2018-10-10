#include <EOSCreator.H>

#include <iostream>

#include <dictionary.H>


using namespace std;


EOS* EOSCreator::create( const string& dictName, const string& eqName ) {


    // Initialize types

    _eosMapType["vanDerWaals"]           = eosType::VdW;
    _eosMapType["Carnahan-Starling"]     = eosType::CS;

    
    // Load model name from dictionary

    dictionary dict(dictName);

    string etype = dict.lookUp<string>( eqName + "/Forces/EOS/type" );


    if( _eosMapType.find(etype) != _eosMapType.end() ) {

	
    	switch( _eosMapType[etype] ) {

	    
    	case eosType::VdW:

    	    return new vanDerWaals(dictName, eqName);

    	    break;



    	case eosType::CS:

    	    return new CarnahanStarling(dictName, eqName);

    	    break;	    

    	}

	
    }

    
    else {


    	// Default
    
    	cout << endl << " [ERROR]  EOS " << etype << " not available" << endl << endl;

    	exit(1);
    
    }

    
    return 0;

}
