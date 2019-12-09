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

#include <writeBasicMesh.H>

#include <periodicBoundaryCorrection.H>



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
	
	vector<string> bdnames = propDict.bracedEntriesNames( "Geometry/boundary" );


	// Polyhedrons for each boundary

	vector< pair<string, Polyhedron> > bdPolyMap;

	for( auto bd : bdnames ) {

	    vector< vector<double> > bbox;

	    cout << "Reading surface " << bd << endl;
	    
	    bdPolyMap.push_back( std::make_pair(  bd, STLToPolyhedron( propDict.lookUp<string>("Geometry/boundary/"+bd+"/file"), bbox ))  ) ;

	}


	// Check for closest boundary

	findClosestBoundary( boundaries, meshPoints, nb, bdPolyMap );
	
    }



    // Periodic correction

    cout << endl << "Applying periodic correction" << endl << endl;

    {

	// Read pairs

	map< pair<string,string>, vector<scalar> > periodicPairs;

	vector<string> bdnames = propDict.bracedEntriesNames( "periodicPairs" );

	for( auto bd : bdnames ) {

	    string other = propDict.lookUp<string>( "periodicPairs/" + bd + "/periodicBoundary" );

	    periodicPairs[ make_pair(bd,other) ] = propDict.lookUp< vector<scalar> >( "periodicPairs/" + bd + "/direction" );

	}


	// Apply correction

	periodicBoundaryCorrection( nb, periodicPairs, meshPoints, boundaries );

    }



    


    // Write final mesh

    cout << endl << "Writing mesh" << endl << endl;

    writeBasicMesh( meshPoints, nb, meshCells, boundaries );



    
    cout << endl << "Finished meshing in " << 0 << " seconds" << endl << endl;

   
    
    return 0;

}
