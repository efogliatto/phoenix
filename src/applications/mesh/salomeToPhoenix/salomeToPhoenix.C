/*

  salomeToPhoenix

  Compute neighbours and periodic boundaries from salome hexahedral mesh

 */


#include <iostream>

#include "readBasicMesh.H"

#include <dictionary.H>

// #include "../meshInclude/latticeMesh_C.H"

#include <latticeModelCreator.H>

// #include "writeLatticeMesh.H"

// #include "computeVirtualNodes.H"



using namespace std;


int main(int argc, char** argv) {



    cout << "                    " << endl;
    cout << "     o-----o-----o  " << endl;
    cout << "     | -   |   - |  " << endl;
    cout << "     |   - | -   |  salomeToPhoenix" << endl;
    cout << "     o<----o---->o  " << endl;
    cout << "     |   - | -   |  Neighbours computing from SALOME mesh" << endl;
    cout << "     | -   |   - |  " << endl;
    cout << "     o-----o-----o  " << endl << endl;



    // Read full mesh

    basicMesh mesh = readBasicMesh();
       

    
    // Read model name. Use D2Q9 as default

    dictionary ldict("properties/latticeProperties");
    
    latticeModelCreator lbm;
	
    latticeModel* lbmodel = lbm.create(ldict.lookUpOrDefault<string>("LBModel","D2Q9"));




    // Compute neighbours

    mesh.nb.resize(mesh.nPoints);

    for(uint i = 0 ; i < mesh.nPoints ; i++)
	mesh.nb[i].resize( lbmodel->q() );


    
    
    



    cout << "Finished mesh transformation" << endl << endl;

   
    
    return 0;

}
