#include <latticeModelCreator.H>
#include <iostream>

using namespace std;

latticeModel* latticeModelCreator::create( const std::string& modelName ) {

    if( modelName.compare("D2Q9") == 0 ) {

	return new D2Q9();

    }

    else {


	if( modelName.compare("D3Q15") == 0 ) {
	    
	    return new D3Q15();

	}

	else {

	    // Default
    
	    cout << endl << "LBModel " << modelName << " does not exist" << endl << endl;

	    exit(1);

	}
    
    }

    
    return 0;

}



latticeModel* latticeModelCreator::create( const uint& d, const uint& q ) {

    if( (d == 2)  &&  (q == 9) ) {

	return new D2Q9();

    }

    else {


	if( (d == 3)  &&  (q == 15) ) {

	    return new D3Q15();

	}

	else {

	    // Default
    
	    cout << endl << "LBModel D" << d << "Q" << q << " does not exist" << endl << endl;

	    exit(1);

	}
    
    }

    
    return 0;

}
