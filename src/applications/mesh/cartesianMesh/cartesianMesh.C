/*

  latticeMeshPartition

  Mesh subdivision for parallel processing

 */


#include <iostream>

#include <latticeModelCreator.H>

#include <dictionary.H>

#include <STLToPolyhedron.H>

#include <STLToNefPolyhedron.H>

#include <createLatticeGrid.H>

#include <createLatticeCells.H>



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
    cout << "                    " << endl;    



    
    // Read model name. Use D2Q9 as default

    dictionary propDict("properties/cartesianMeshProperties");
    
    latticeModelCreator lbm;
	
    latticeModel* lbmodel = lbm.create(propDict.lookUpOrDefault<string>("LBModel","D2Q9"));



    // Load Main Geometry

    vector< vector<double> > bbox = {{0,0,0}, {0,0,0}};
    
    Nef_polyhedron P = STLToNefPolyhedron(  propDict.lookUp<string>("Geometry/name"), bbox );



    
    // Build base mesh

    uint nx( bbox[1][0] - bbox[0][0] );

    uint ny( bbox[1][1] - bbox[0][1] );

    uint nz( bbox[1][2] - bbox[0][2] );


    cout << endl << "Building base grid with " << nx*ny*nz << " points" << endl;
	
    vector< vector<uint> > basePoints;

    createLatticeGrid( basePoints, lbmodel, nx, ny, nz );



    cout << endl << "Building base cells with " << (nx-1)*(ny-1)*(nz-1) << " cells" << endl;
	
    vector< vector<uint> > baseCells;

    createLatticeCells( baseCells, lbmodel, nx, ny, nz );    


    

    // Check which cell centers are inside Polyhedron

    vector<uint> isInside;

    
    
    



    cout << endl << "Finished meshing" << endl << endl;

   
    
    return 0;

}
