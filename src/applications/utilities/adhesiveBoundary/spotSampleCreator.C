#include <spotSampleCreator.H>

#include <dictionary.H>

using namespace std;


/** Model creation */
    
unique_ptr<spotSample> spotSampleCreator::create( const string& sptype, const string& bdname ) {


    // Initialize types

    _spotMapType["uniform"] = spotType::uniform;
    
    _spotMapType["random"]  = spotType::random;        



    
    // // Load model name from dictionary

    // dictionary dict("properties/macroProperties");

    // string modelName = dict.lookUpOrDefault<string>( entry + "/relaxModel", "uniform" );


    

    // Assign model
    
    if( _spotMapType.find(sptype) != _spotMapType.end() ) {

	
    	switch( _spotMapType[sptype] ) {
	    

    	case spotType::uniform:

	    return unique_ptr<spotSample>( new uniformSpots(bdname) );

    	    break;


    	case spotType::random:

	    return unique_ptr<spotSample>( new randomSpots(bdname) );

    	    break;


    	}

    }


    else {
    
    	cout << endl << " [ERROR]  Spot types " << sptype << " not available" << endl << endl;

    	exit(1);

    }    
    
    

    return 0;

}
