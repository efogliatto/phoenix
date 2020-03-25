/*

  latticeBox3D

  Lattice generation in a 3D Box

 */


#include <iostream>

#include <dictionary.H>

#include <latticeModelCreator.H>

#include "writeBasicMesh.H"

#include "genericBoundary3D.H"

#include "periodicXY3D.H"

#include "periodicXYZ3D.H"

#include "periodicZ3D.H"




using namespace std;


int main(int argc, char** argv) {



    cout << "                    " << endl;
    cout << "     o-----o-----o  " << endl;
    cout << "     | -   |   - |  " << endl;
    cout << "     |   - | -   |   latticeBox3D" << endl;
    cout << "     o<----o---->o  " << endl;
    cout << "     |   - | -   |  3D lattice box" << endl;
    cout << "     | -   |   - |  " << endl;
    cout << "     o-----o-----o  " << endl << endl;



    // Read full mesh

    basicMesh mesh;
    

    
    
    // Grid size

    dictionary ldict("properties/latticeProperties");

    uint nx( (uint)ldict.lookUp<scalar>("Nx") );

    uint ny( (uint)ldict.lookUp<scalar>("Ny") );

    uint nz( (uint)ldict.lookUp<scalar>("Nz") );    

    mesh.nPoints = nx*ny*nz;

   

    
    // Read model name. Use D3Q15 as default
    
    latticeModelCreator lbm;
	
    latticeModel* lbmodel = lbm.create(ldict.lookUpOrDefault<string>("LBModel","D3Q15"));

    mesh.D = lbmodel->d();

    mesh.Q = lbmodel->q();    
    
    


    
    // ******************************************************************** //
    //                         Points inside geometry                       //
    // ******************************************************************** //

    cout << "Adding points to lattice" << endl << endl;
    
    mesh.points.resize( mesh.nPoints );

    for(uint i = 0 ; i < mesh.nPoints ; i++)
    	mesh.points[i].resize(3);


    for( uint k = 0 ; k < nz ; k++) {
    
	for( uint j = 0 ; j < ny ; j++) {
	
	    for( uint i = 0 ; i < nx ; i++) {

		uint idx = i + j*nx + k*nx*ny;
			    
		mesh.points[idx][0] = i;
		
		mesh.points[idx][1] = j;

		mesh.points[idx][2] = k;
		
	    }
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
    // p with index pointId is at most at pointId + x + (y*Nx) + (z*Nx*Ny).
    // For boundary nodes, check neighbouring using distance to point
    

    for( uint k = 0 ; k < nz ; k++ ) {

    	for( uint j = 0 ; j < ny ; j++) {
	
    	    for( uint i = 0 ; i < nx ; i++) {


		uint idx = i + j*nx + k*nx*ny;
		
		
		// Boundary nodes
		
    		if( (i == 0)  ||  (i == nx-1)   ||   (j == 0)  ||  (j == ny-1)   ||   (k == 0)  ||  (k == nz-1) ) {	       		    

    		    for( uint vid = 0 ; vid < mesh.Q ; vid++ ) {

    			int newId = idx   +   vel[ reverse[vid] ][0]   +   vel[ reverse[vid] ][1] * nx   +   vel[ reverse[vid] ][2] * nx * ny;

    			if( newId >= 0   &&   newId <= (int)(nx*ny*nz-1) ) {

    			    if (      (  abs( mesh.points[idx][0] - mesh.points[newId][0] ) <= 1  )
    			         &&   (  abs( mesh.points[idx][1] - mesh.points[newId][1] ) <= 1  )
    			         &&   (  abs( mesh.points[idx][2] - mesh.points[newId][2] ) <= 1  )  ) {
	    
    				mesh.nb[idx][vid] = newId;

    			    }

    			}

    		    }

    		}


		
		// Inner nodes
		
		else {

		    for( uint vid = 0 ; vid < mesh.Q ; vid++ ) {

			mesh.nb[idx][vid] = idx   +   vel[ reverse[vid] ][0]   +   vel[ reverse[vid] ][1] * nx    +    vel[ reverse[vid] ][2] * nx * ny;

		    }

		}

		
	    
    	    }

    	}

    }    
    








    // ******************************************************************** //
    //                             Boundary                                 //
    // ******************************************************************** //

    cout << "Computing boundary nodes" << endl << endl;

    
    // Boundary type
    
    string bdt = ldict.lookUpOrDefault<string>( "boundaryType", "generic");

    
    // Assign points on boundary based on bdType

    // Generic
    if( bdt == "generic") {

    	genericBoundary3D( mesh, nx, ny, nz );

	mesh.bd.bdNames.push_back( "X0" );
	
	mesh.bd.bdNames.push_back( "X1" );
	
	mesh.bd.bdNames.push_back( "Y0" );
	
	mesh.bd.bdNames.push_back( "Y1" );

	mesh.bd.bdNames.push_back( "Z0" );
	
	mesh.bd.bdNames.push_back( "Z1" );	

    }


    // periodicXY
    else {

    	if( bdt == "periodicXY" ) {

    	    periodicXY3D( mesh, nx, ny, nz );

	    mesh.bd.bdNames.push_back( "Z0" );
	
	    mesh.bd.bdNames.push_back( "Z1" );	    

    	}


    	// periodicZ
    	else {

    	    if( bdt == "periodicZ" ) {

    		periodicZ3D( mesh, nx, ny, nz );

		mesh.bd.bdNames.push_back( "X0" );
	
		mesh.bd.bdNames.push_back( "X1" );
	
		mesh.bd.bdNames.push_back( "Y0" );
	
		mesh.bd.bdNames.push_back( "Y1" );

    	    }


    	    // Fully periodic
    	    else {

    		if( bdt == "periodicXYZ" ) {

    		    periodicXYZ3D( mesh, nx, ny, nz );

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
    
    mesh.vtkCells.resize( (nx-1)*(ny-1)*(nz-1) );

    for( uint i = 0 ; i < mesh.vtkCells.size() ; i++ )
    	mesh.vtkCells[i].resize( 8 );

    
    mesh.ncells = 0;
    
    mesh.cellType = 8;

    for( uint k = 0 ; k < nz-1 ; k++ ) {
    
	for( uint j = 0 ; j < ny-1 ; j++ ) {
	
	    for( uint i = 0 ; i < (nx-1) ; i++ ) {

		mesh.vtkCells[mesh.ncells][0] = i   +   j*nx             +   k*nx*ny;

		mesh.vtkCells[mesh.ncells][1] = i   +   j*nx + 1         +   k*nx*ny;

		mesh.vtkCells[mesh.ncells][2] = i   +   j*nx + nx        +   k*nx*ny;

		mesh.vtkCells[mesh.ncells][3] = i   +   j*nx + nx + 1    +   k*nx*ny;

		mesh.vtkCells[mesh.ncells][4] = i   +   j*nx             +   (k+1)*nx*ny;

		mesh.vtkCells[mesh.ncells][5] = i   +   j*nx + 1         +   (k+1)*nx*ny;

		mesh.vtkCells[mesh.ncells][6] = i   +   j*nx + nx        +   (k+1)*nx*ny;

		mesh.vtkCells[mesh.ncells][7] = i   +   j*nx + nx + 1    +   (k+1)*nx*ny;
	
		mesh.ncells++;

	    }

	}

    }
    
    
    



    // Write mesh

    cout << "Writing mesh" << endl << endl;

    writeBasicMesh( mesh );

    

    



    cout << "Finished meshing" << endl << endl;
    
   
    
    return 0;

}
