#include <cavityModelCreator.H>

#include <iostream>

using namespace std;


/** Model creation */
    
unique_ptr<cavityModel> cavityModelCreator::create( const string& bdname ) {
   
    
    // Initialize types

    _cavityMapType["flat"]   = cavityType::flat;
    
    _cavityMapType["linear"] = cavityType::linear;        

    

    // Read type for this boundary

    dictionary dict("properties/adhesiveProperties");

    string cavtype = dict.lookUpOrDefault<string>("Boundaries/" + bdname + "/CavityModel/type", "flat");

    

    // Assign model
    
    if( _cavityMapType.find(cavtype) != _cavityMapType.end() ) {

	
    	switch( _cavityMapType[cavtype] ) {
	    

    	case cavityType::flat:

	    return unique_ptr<cavityModel>( new flatCavity(bdname) );

    	    break;


    	case cavityType::linear:

	    return unique_ptr<cavityModel>( new linearCavity(bdname) );

    	    break;


    	}

    }


    else {
    
    	cout << endl << " [ERROR]  Spot types " << cavtype << " not available" << endl << endl;

    	exit(1);

    }    
    
    

    return 0;

}
