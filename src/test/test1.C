#include <iostream>
#include <latticeModelCreator.H>


using namespace std;

int main() {

    latticeModelCreator lbcreator;
    
    latticeModel* lbm = lbcreator.create("D2Q9");

    for(uint i = 0 ; i < lbm->q() ; i++)
	cout << lbm->lvel()[i][0] << "  " << lbm->lvel()[i][1] << "  " << lbm->lvel()[i][2] << "  " << endl;
    
    return 0;

}
