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

#include <cellsInsidePolyhedron.H>

#include <updatePointsAndCells.H>

#include <computeNeighboursFromCells.H>

#include <findClosestBoundary.H>



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

    vector< vector<double> > bbox;
    
    Polyhedron P = STLToPolyhedron(  propDict.lookUp<string>("Geometry/name"), bbox );




    // Build mesh

    vector< vector<uint> > meshCells;

    vector< vector<uint> > meshPoints;    
    

    {
    
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

	cout << endl << "Finding base cell centers inside main geometry" << endl;

	vector<bool> isInside;

	cellsInsidePolyhedron( isInside, P, basePoints, baseCells );
    


	// Clear unused points and cells

	cout << endl << "Reducing mesh:" << endl;

	updatePointsAndCells( basePoints, baseCells, meshPoints, meshCells, isInside );

    
    }


    cout << "Final mesh with " << meshPoints.size() << " points and " << meshCells.size() << " cells" << endl;



    // Neighbours

    cout << endl << "Computing neighbours using cell information" << endl;

    vector< vector<int> > nb;

    computeNeighboursFromCells( nb, meshPoints, meshCells, lbmodel );


    
    // Assign boundaries

    cout << endl << "Boundary assignment" << endl << endl;

    map< string, vector<uint> > boundaries;

    {

	// Boundaries names
	
	vector<string> bdnames = propDict.bracedEntry( "Geometry/boundary" );


	// Polyhedrons for each boundary

	vector< pair<string, Polyhedron> > bdPolyMap;

	for( auto bd : bdnames ) {

	    vector< vector<double> > bbox;
	    
	    bdPolyMap.push_back( std::make_pair(  bd, STLToPolyhedron( bd, bbox ))  ) ;

	}


	// Check for closest boundary

	findClosestBoundary( boundaries, meshPoints, nb, bdPolyMap );
	
    }



    cout << endl << "Finished meshing" << endl << endl;

   
    
    return 0;

}
