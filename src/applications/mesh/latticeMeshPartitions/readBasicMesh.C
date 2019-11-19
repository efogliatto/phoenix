#include "readBasicMesh.H"

#include <iostream>

#include <fstream>

#include <dictionary.H>

#include <latticeModelCreator.H>

using namespace std;



basicMesh readBasicMesh() {

    basicMesh mesh;

    uint i,j;

    uint status;

    

    // ******************************************************************** //
    //                         Points inside geometry                       //
    // ******************************************************************** //

    cout << "Reading Mesh points" << endl << endl;


    // Open file and load content

    ifstream inFile;

    inFile.open( "lattice/points" );

    if( inFile.is_open() == false ){
	
    	cout << " [ERROR]  Unable to find file lattice/points" << endl;
	
    	exit(1);
	
    }

    
    
    // Number of points
    
    int iaux;

    inFile >> iaux;

    mesh.nPoints = (uint)iaux;

    mesh.points.resize(mesh.nPoints);

    for(uint i = 0 ; i < mesh.nPoints ; i++)
	mesh.points[i].resize(3);
    
    
    // Read Mesh points

    for( i = 0 ; i < mesh.nPoints ; i++ ) {

    	inFile >> mesh.points[i][0];
    	inFile >> mesh.points[i][1];
    	inFile >> mesh.points[i][2];

    }

    inFile.close();


    

    // ******************************************************************** //
    //                             Neighbours                               //
    // ******************************************************************** //

    cout << "Reading Neighbour indices" << endl << endl;


    // Open file and load content

    inFile.open( "lattice/neighbours" );

    if( inFile.is_open() == false ){
	
    	cout << " [ERROR]  Unable to find file lattice/neighbours" << endl;
	
    	exit(1);
	
    }



    
    // Read model name. Use D2Q9 as default

    dictionary dict("properties/latticeProperties");
    
    latticeModelCreator lbm;
	
    latticeModel* lbmodel = lbm.create(dict.lookUpOrDefault<string>("LBModel","D2Q9"));

    mesh.Q = lbmodel->q();




    
    // Read neighbours

    mesh.nb.resize(mesh.nPoints);

    for(uint i = 0 ; i < mesh.nPoints ; i++)
	mesh.nb[i].resize(mesh.Q);    
    
    for( uint i = 0 ; i < mesh.nPoints ; i++ ) {

    	for( uint j = 0 ; j < mesh.Q ; j++ ) {

    	    inFile >> mesh.nb[i][j];

    	}

    }

    inFile.close();

    


    
    // ******************************************************************** //
    //                             VTK Cells                                //
    // ******************************************************************** //

    cout << "Reading vtkCells" << endl << endl;


    // Open file and load content

    inFile.open( "lattice/vtkCells" );

    if( inFile.is_open() == false ){
	
    	cout << " [ERROR]  Unable to find file lattice/vtkCells" << endl;
	
    	exit(1);
	
    }

    
    // Number of cells and points per cell
    
    inFile >> iaux;

    mesh.ncells = (uint)iaux;
    
    inFile >> iaux;
    
    mesh.cellType = (uint)iaux;

    
    
    // Read cells

    mesh.vtkCells.resize(mesh.ncells);

    for(uint i = 0 ; i < mesh.ncells ; i++)
	mesh.vtkCells[i].resize(mesh.cellType);
    
    
    for( uint i = 0 ; i < mesh.ncells ; i++ ) {

    	for( uint j = 0 ; j < mesh.cellType ; j++ ) {

    	    inFile >> mesh.vtkCells[i][j];

    	}

    }

    inFile.close();






    // ******************************************************************** //
    //                        Reading boundary nodes                        //
    // ******************************************************************** //

    cout << "Reading boundary nodes" << endl << endl;


    // Open file and load content

    inFile.open( "lattice/boundary" );

    if( inFile.is_open() == false ){
	
    	cout << " [ERROR]  Unable to find file lattice/boundary" << endl;
	
    	exit(1);
	
    }    

   
    // Number of boundary types
    
    inFile >> iaux;

    mesh.bd.nbd = (uint)iaux;

    
    // Total number of elements per boundary type

    mesh.bd.nbdelem.resize( mesh.bd.nbd );

    mesh.bd.bdNames.resize( mesh.bd.nbd );

    
    // Elements in boundary

    mesh.bd.bdPoints.resize( mesh.bd.nbd );
    
   
    // Read boundary
    
    for( uint i = 0 ; i < mesh.bd.nbd ; i++ ) {

	
    	// Boundary name
	
    	inFile >> mesh.bd.bdNames[i];

	
    	// Elements in boundary

	inFile >> iaux;

	mesh.bd.nbdelem[i] = (uint)iaux;

	
    	// Resize bdPoints
	
    	mesh.bd.bdPoints[i].resize( mesh.bd.nbdelem[i] );
	
    	for( uint j = 0 ; j < mesh.bd.nbdelem[i] ; j++ ) {

    	    inFile >> iaux;
	    
    	    mesh.bd.bdPoints[i][j] = (uint)iaux;

    	}

    }
    

    inFile.close();

    if (status) {}
    
    return mesh;
    
}
