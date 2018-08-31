#include <LBModelCreator.h>

using namespace std;

basicLBModel* LBModelCreator::create(const std::string& modelName) {

    if( modelName.compare("D2Q9") == 0 ) return new D2Q9();
    if( modelName.compare("D2Q5") == 0 ) return new D2Q5();
    if( modelName.compare("D2Q4") == 0 ) return new D2Q4();
    if( modelName.compare("D3Q7") == 0 ) return new D3Q7();


    // Default 
    cout << endl << "LBModel " << modelName << " does not exist" << endl << endl;
    exit(1);

    return 0;

}
