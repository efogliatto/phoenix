#include <iostream>
#include <latticeModelCreator.H>


using namespace std;

int main() {

    latticeModelCreator lbcreator;
    
    latticeModel* lbm = lbcreator.create("D3Q15");

    const scalarMatrix& M = lbm->MRTMatrix();

    const scalarMatrix& invM = lbm->MRTInvMatrix(); 

    for(uint i = 0 ; i < lbm->q() ; i++) {

	for(uint j = 0 ; j < lbm->q() ; j++) {

	    scalar a = 0;

	    for(uint k = 0 ; k < lbm->q() ; k++) {

		a += M[i][k] * invM[k][j];

	    }

	    cout << a << "  ";

	}

	cout << endl;

    }
    
    return 0;

}
