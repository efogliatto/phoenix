/*

  latticeBox2D

  lattice generation in a 2D Box

*/

#include <iostream>

#include "basicMesh.H"

#include <dictionary.H>

#include <latticeModelCreator.H>



// void periodicX( basicMesh* mesh, uint nx, uint ny );
// void periodicY( basicMesh* mesh, uint nx, uint ny );
// void periodicXY( basicMesh* mesh, uint nx, uint ny );
// void genericBoundary( basicMesh* mesh, uint nx, uint ny );



using namespace std;


int main(int argc, char** argv) {

    
    cout << "\n\n   LATTICE MESH GENERATION\n\n";


    
    // Basic lattice properties

    basicMesh mesh;

    

    // Lattice size

    uint nx = 0,
    	ny = 0;

    dictionary latprop("properties/latticeProperties");

    nx = (uint)latprop.lookUp<scalar>("Nx");

    ny = (uint)latprop.lookUp<scalar>("Ny");

    mesh.nPoints = nx*ny;




    
    // Lattice model

    latticeModelCreator lbcreator;
    
    latticeModel* lbm = lbcreator.create( latprop.lookUp<string>("LBModel") );
    

    

    

    // ******************************************************************** //
    //                         Points inside geometry                       //
    // ******************************************************************** //

    
    cout << "Adding points to lattice\n\n";

    mesh.points.resize( mesh.nPoints );

    for(uint i = 0 ; i < mesh.nPoints ; i++)
	mesh.points[i].resize(3);
    

        
    // Move over indices

    for( uint j = 0 ; j < ny ; j++) {
	
    	for( uint i = 0 ; i < nx ; i++) {
			    
    	    mesh.points[i+j*nx][0] = i;
    	    mesh.points[i+j*nx][1] = j;
		
    	}
    }




    
    

    // ******************************************************************** //
    //                             Neighbours                               //
    // ******************************************************************** //

    
    cout << "Computing neighbour indices\n\n";

    
    // Lattice velocities
    
    const vector< vector<int> >& velocities = lbm->lvel();

    const vector<uint>& rev = lbm->reverse();

    
    
    // Create and resize neighbour matrix
    
    mesh.Q = lbm->q();

    mesh.nb.resize( mesh.nPoints );

    for(uint i = 0 ; i < mesh.nPoints ; i++)
	mesh.nb[i].resize(mesh.Q);    

    for(uint i = 0 ; i < mesh.nPoints ; i++) {
	
	for(uint j = 0 ; j < mesh.Q ; j++) {
	    
	    mesh.nb[i][j] = -1;

	}

    }

    

    

    // Check for neighbouring
    // There is no need to iterate over all points to look for a neighbour. Given a lattice velocity vector (x,y,z), the neighbour of a point
    // p with index pointId is at most at pointId + x + (y*Nx).

    
    // Internal points first

    int pointId;
    
    for( uint j = 1 ; j < (ny-1) ; j++) {
	
    	for( uint i = 1 ; i < (nx-1) ; i++) {

    	    pointId = i+j*nx;


    	    // Iterate on velocities
    	    for( int velId = 0 ; velId < (int)mesh.Q ; velId++ ) {

    		mesh.nb[pointId][velId] = pointId   +   velocities[ rev[velId] ][0]   +   velocities[ rev[velId] ][1] * nx;

    	    }
	    
    	}

    }

    

    // For boundary nodes, check neighbouring using distance to point

    for( uint j = 0 ; j < ny ; j+=(ny-1)) {
	
    	for( uint i = 0 ; i < nx ; i++) {

    	    pointId = i+j*nx;

    	    // Iterate on velocities
    	    for( int velId = 0 ; velId < (int)mesh.Q ; velId++ ) {

    		int newId = pointId   +   velocities[ rev[velId] ][0]   +   velocities[ rev[velId] ][1] * nx;

    		if( newId >= 0   &&   newId <= (int)(nx*ny-1) ) {

    		    if ( (  abs( mesh.points[pointId][0] - mesh.points[newId][0] ) <= 1  )   &&   (  abs( mesh.points[pointId][1] - mesh.points[newId][1] ) <= 1  )  ) {
	    
    			mesh.nb[pointId][velId] = newId;

    		    }

    		}

    	    }
	    
    	}

    }


    for( uint j = 1 ; j < ny-1 ; j++ ) {
	
    	for( uint i = 0 ; i < nx ; i+=(nx-1)) {

    	    pointId = i+j*nx;

    	    // Iterate on velocities
    	    for( int velId = 0 ; velId < (int)mesh.Q ; velId++ ) {

    		int newId = pointId   +   velocities[ rev[velId] ][0]   +   velocities[ rev[velId] ][1] * nx;

    		if( newId >= 0   &&   newId <= (int)(nx*ny-1) ) {

    		    if ( (  abs( mesh.points[pointId][0] - mesh.points[newId][0] ) <= 1  )   &&   (  abs( mesh.points[pointId][1] - mesh.points[newId][1] ) <= 1  )  ) {
	    
    			mesh.nb[pointId][velId] = newId;

    		    }

    		}

    	    }
	    
    	}

    }
    


    	
	    
    
    // ******************************************************************** //
    //                             Boundary                                 //
    // ******************************************************************** //

    cout << "Computing boundary nodes\n\n";

    mesh.bd.bdNames.push_back("X0");
    mesh.bd.bdNames.push_back("X1");
    mesh.bd.bdNames.push_back("Y0");
    mesh.bd.bdNames.push_back("Y1");

    
    // // Boundary type
    
    // char* bdt;

    // status = lookUpStringEntry("properties/latticeProperties","boundaryType", &bdt, "generic");

    
    // // Assign points on boundary based on bdType

    // // Generic
    // if( strcmp(bdt,"generic") == 0) {

    // 	genericBoundary( &mesh, nx, ny );

    // }


    // // periodicX
    // else {

    // 	if( strcmp(bdt,"periodicX")  == 0 ) {

    // 	    periodicX( &mesh, nx, ny );

    // 	}


    // 	// periodicY
    // 	else {

    // 	    if( strcmp(bdt,"periodicY")  == 0 ) {

    // 		periodicY( &mesh, nx, ny );

    // 	    }


    // 	    // periodicXY
    // 	    else {

    // 		if( strcmp(bdt,"periodicXY")  == 0 ) {

    // 		    periodicXY( &mesh, nx, ny );

    // 		}

    // 		else {

    // 		    printf("[ERROR]   Unrecognized boundary type %s\n\n", bdt);
    // 		    exit(0);

    // 		}
	
    // 	    }
	
    // 	}

	
    // }
    




	    

    


    // ******************************************************************** //
    //                             VTK Cells                                //
    // ******************************************************************** //

    cout << "Creating VTK Cells\n\n";
    
    mesh.vtkCells.resize( (nx-1)*(ny-1) );

    for(uint i = 0 ; i < mesh.vtkCells.size() ; i++)
	mesh.vtkCells[i].resize(4);        

    
    mesh.ncells = 0;
    mesh.cellType = 4;
    
    for( uint j = 0 ; j < ny-1 ; j++ ) {
	
    	for( uint i = 0 ; i < (nx-1) ; i++ ) {

    	    mesh.vtkCells[mesh.ncells][0] = i + j*nx;
    	    mesh.vtkCells[mesh.ncells][1] = i + j*nx + 1;
    	    mesh.vtkCells[mesh.ncells][2] = i + j*nx + nx;
    	    mesh.vtkCells[mesh.ncells][3] = i + j*nx + nx + 1;
	
    	    mesh.ncells++;

    	}

    }
    

    



    


	    


    // // ******************************************************************** //
    // //                             Writing                                  //
    // // ******************************************************************** //
    
    // printf("Writting Mesh\n\n");

    // writeBasicMesh( &mesh );


    


    

    return 0;

}
