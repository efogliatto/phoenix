/*

  latticeMeshPartition

  Mesh subdivision for parallel processing

 */


#include <iostream>

#include <latticeModelCreator.H>

#include <dictionary.H>




using namespace std;


int main(int argc, char** argv) {



    cout << "                    " << endl;
    cout << "     o-----o-----o  " << endl;
    cout << "     | -   |   - |  " << endl;
    cout << "     |   - | -   |        cartesianMesh" << endl;
    cout << "     o<----o---->o  " << endl;
    cout << "     |   - | -   |  Cartesian Mesh Generator" << endl;
    cout << "     | -   |   - |  " << endl;
    cout << "     o-----o-----o  " << endl << endl;



    
    // Read model name. Use D2Q9 as default

    dictionary ldict("properties/latticeProperties");
    
    latticeModelCreator lbm;
	
    latticeModel* lbmodel = lbm.create(ldict.lookUpOrDefault<string>("LBModel","D2Q9"));

    



    cout << "Finished meshing" << endl << endl;

   
    
    return 0;

}
