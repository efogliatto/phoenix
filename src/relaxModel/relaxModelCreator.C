#include <relaxModelCreator.H>

#include <dictionary.H>


using namespace std;


/** Model creation */

unique_ptr<relaxModel> relaxModelCreator::create( const std::string& entry ) {


    
    // Initialize types

    _relaxMapType["uniform"]            = relaxType::utau;
    
    _relaxMapType["rhoPieceWise"]       = relaxType::rhopw;

    _relaxMapType["rhoPieceWiseLinear"] = relaxType::rhopwl;    



    
    // Load model name from dictionary

    dictionary dict("properties/macroProperties");

    string modelName = dict.lookUpOrDefault<string>( entry + "/relaxModel", "uniform" );


    

    // Assign model
    
    if( _relaxMapType.find(modelName) != _relaxMapType.end() ) {

	switch( _relaxMapType[modelName] ) {

	case relaxType::utau:

	    return std::make_unique<uniformTau>(entry);

	    break;


	case relaxType::rhopw:

	    return std::make_unique<rhoPieceWise>(entry);

	    break;


	case relaxType::rhopwl:

	    return std::make_unique<rhoPieceWiseLinear>(entry);

	    break;	    

	}

    }


    else {
    
	cout << endl << " [ERROR]  Relaxation model " << modelName << " not available" << endl << endl;

	exit(1);

    }    
    
    

    return 0;

}
