#include <EOSCreator.H>

#include <iostream>

#include <dictionary.H>


using namespace std;


EOS* EOSCreator::create( const std::string& fname ) {

    
    // Load model name from dictionary

    dictionary dict(fname);

    string model = dict.lookUp<string>("EOS/model");

    
    if( model.compare("vanDerWaals") == 0 ) {

    	return new vanDerWaals(fname);

    }

    else {


    	// Default
    
    	cout << endl << " [ERROR]  EOS " << model << " not available" << endl << endl;

    	exit(1);
    
    }

    
    return 0;

}
