#include <spotRadiusCreator.H>

using namespace std;


/** Model creation */
    
unique_ptr<spotRadius> spotRadiusCreator::create( const string& rtype ) {

    
    // Initialize types

    _radiusMapType["fixed"]  = radiusType::fixed;
    
    _radiusMapType["normal"] = radiusType::normal;        


   

    // Assign model
    
    if( _radiusMapType.find(rtype) != _radiusMapType.end() ) {

	
    	switch( _radiusMapType[rtype] ) {
	    

    	case radiusType::fixed:

	    return unique_ptr<spotRadius>( new fixedRadiusSpot() );

    	    break;


    	case radiusType::normal:

	    return unique_ptr<spotRadius>( new normalDistRadiusSpot() );

    	    break;


    	}

    }


    else {
    
    	cout << endl << " [ERROR]  Spot radius type " << rtype << " not available" << endl << endl;

    	exit(1);

    }    
    
    

    return 0;

}
