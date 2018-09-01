#include <latticeModelCreator.H>
#include <iostream>

using namespace std;

latticeModel* latticeModelCreator::create( const std::string& modelName ) {

    if( modelName.compare("D2Q9") == 0 ) {

	return new D2Q9();

    }

    else {


	// Default
    
	cout << endl << "LBModel " << modelName << " does not exist" << endl << endl;

	exit(1);
    
    }

    
    return 0;

}
