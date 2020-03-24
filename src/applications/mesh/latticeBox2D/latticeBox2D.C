/*

  latticeBox2D

  Lattice generation in a 2D Box

 */


#include <iostream>

#include <dictionary.H>

#include <latticeModelCreator.H>

#include "writeBasicMesh.H"

#include "genericBoundary.H"

#include "periodicX.H"

#include "periodicY.H"

#include "periodicXY.H"




using namespace std;


int main(int argc, char** argv) {



    cout << "                    " << endl;
    cout << "     o-----o-----o  " << endl;
    cout << "     | -   |   - |  " << endl;
    cout << "     |   - | -   |   latticeBox2D" << endl;
    cout << "     o<----o---->o  " << endl;
    cout << "     |   - | -   |  2D lattice box" << endl;
    cout << "     | -   |   - |  " << endl;
    cout << "     o-----o-----o  " << endl << endl;



    // Read full mesh

    basicMesh mesh;
    

    
    
    // Grid size

    dictionary ldict("properties/latticeProperties");

    uint nx( (uint)ldict.lookUp<scalar>("Nx") );

    uint ny( (uint)ldict.lookUp<scalar>("Ny") );

    mesh.nPoints = nx*ny;

   

    
    // Read model name. Use D2Q9 as default
    
    latticeModelCreator lbm;
	
    latticeModel* lbmodel = lbm.create(ldict.lookUpOrDefault<string>("LBModel","D2Q9"));

    mesh.D = lbmodel->d();

    mesh.Q = lbmodel->q();    
    
    


    
    // ******************************************************************** //
    //                         Points inside geometry                       //
    // ******************************************************************** //

    cout << "Adding points to lattice" << endl << endl;
    
    mesh.points.resize( mesh.nPoints );

    for(uint i = 0 ; i < mesh.nPoints ; i++)
    	mesh.points[i].resize(3);

    
    for( uint j = 0 ; j < ny ; j++ ) {
	
    	for( uint i = 0 ; i < nx ; i++ ) {
			    
    	    mesh.points[i+j*nx][0] = i;
    	    mesh.points[i+j*nx][1] = j;
		
    	}
    }





    // ******************************************************************** //
    //                             Neighbours                               //
    // ******************************************************************** //

    cout << "Computing neighbour indices" << endl << endl;


    // Lattice constants
    
    const vector< vector<int> > vel = lbmodel->lvel();

    vector<uint> reverse = lbmodel->reverse();


    // Neighbour array

    mesh.nb.resize( mesh.nPoints );

    for( uint i = 0 ; i < mesh.nPoints ; i++ ) {

    	mesh.nb[i].resize( mesh.Q );

    	std::fill( mesh.nb[i].begin(), mesh.nb[i].end(), -1 );

    }




    // Check for neighbouring
    // There is no need to iterate over all points to look for a neighbour. Given a lattice velocity vector (x,y,z), the neighbour of a point
    // p with index pointId is at most at pointId + x + (y*Nx).

    
    // Internal points first
    
    for( uint j = 1 ; j < (ny-1) ; j++) {
	
    	for( uint i = 1 ; i < (nx-1) ; i++) {

    	    uint pointId = i+j*nx;


    	    // Iterate on velocities
	    
    	    for( uint velId = 0 ; velId < mesh.Q ; velId++ )
    		mesh.nb[pointId][velId] = pointId   +   vel[ reverse[velId] ][0]   +   vel[ reverse[velId] ][1] * nx;

	    
    	}

    }

    

    // For boundary nodes, check neighbouring using distance to point

    for( uint j = 0 ; j < ny ; j+=(ny-1)) {
	
    	for( uint i = 0 ; i < nx ; i++) {

    	    uint pointId = i+j*nx;

	    
    	    // Iterate on velocities
	    
    	    for( uint velId = 0 ; velId < mesh.Q ; velId++ ) {

    		int newId = pointId   +   vel[ reverse[velId] ][0]   +   vel[ reverse[velId] ][1] * nx;

    		if( newId >= 0   &&   newId <= nx*ny-1 ) {

    		    if ( (  abs( mesh.points[pointId][0] - mesh.points[newId][0] ) <= 1  )   &&   (  abs( mesh.points[pointId][1] - mesh.points[newId][1] ) <= 1  )  ) {
	    
    			mesh.nb[pointId][velId] = newId;

    		    }

    		}

    	    }
	    
    	}

    }
    


    for( uint j = 1 ; j < ny-1 ; j++ ) {
	
    	for( uint i = 0 ; i < nx ; i+=(nx-1)) {

    	    uint pointId = i+j*nx;

	    
    	    // Iterate on velocities
	    
    	    for( uint velId = 0 ; velId < mesh.Q ; velId++ ) {

    		int newId = pointId   +   vel[ reverse[velId] ][0]   +   vel[ reverse[velId] ][1] * nx;

    		if( newId >= 0   &&   newId <= nx*ny-1 ) {

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

    cout << "Computing boundary nodes" << endl << endl;

    mesh.bd.bdNames.push_back( "X0" );
    mesh.bd.bdNames.push_back( "X1" );
    mesh.bd.bdNames.push_back( "Y0" );
    mesh.bd.bdNames.push_back( "Y1" );

    
    // Boundary type
    
    string bdt = ldict.lookUpOrDefault<string>( "boundaryType", "generic");

    
    // Assign points on boundary based on bdType

    // Generic
    if( bdt == "generic") {

    	genericBoundary( mesh, nx, ny );

    }


    // periodicX
    else {

    	if( bdt == "periodicX" ) {

    	    periodicX( mesh, nx, ny );

    	}


    	// periodicY
    	else {

    	    if( bdt == "periodicY" ) {

    		periodicY( mesh, nx, ny );

    	    }


    	    // periodicXY
    	    else {

    		if( bdt == "periodicXY" ) {

    		    periodicXY( mesh, nx, ny );

    		}

    		else {

    		    cout << " [ERROR]   Unrecognized boundary type " << bdt << endl << endl;
		    
    		    exit(1);

    		}
	
    	    }
	
    	}

	
    }








    



    // ******************************************************************** //
    //                             VTK Cells                                //
    // ******************************************************************** //

    cout << "Creating VTK Cells" << endl << endl;
    
    mesh.vtkCells.resize( (nx-1)*(ny-1) );

    for( uint i = 0 ; i < mesh.vtkCells.size() ; i++ )
    	mesh.vtkCells[i].resize( 4 );

    
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

    
    



    // Write mesh

    cout << "Writing mesh" << endl << endl;

    writeBasicMesh( mesh );

    

    



    cout << "Finished meshing" << endl << endl;
    
   
    
    return 0;

}
