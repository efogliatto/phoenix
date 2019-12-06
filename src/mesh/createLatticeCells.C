#include <createLatticeCells.H>

using namespace std;


void createLatticeCells( vector< vector<uint> >& cells, const latticeModel* lbmodel, const uint nx, const uint ny, const uint nz ) {

    
    // Create grid according to lattice model


    // 2D
    
    if( lbmodel->d() == 2 ) {

	cells.resize( (nx-1)*(ny-1) );

	for(uint i = 0 ; i < (nx-1)*(ny-1) ; i++)
	    cells[i].resize(4,0);


	uint ncells(0);
	
	for( uint j = 0 ; j < (ny-1) ; j++ ) {
	
	    for( uint i = 0 ; i < (nx-1) ; i++ ) {

		cells[ncells][0] = i + j*nx;
		cells[ncells][1] = i + j*nx + 1;
		cells[ncells][2] = i + j*nx + nx;
		cells[ncells][3] = i + j*nx + nx + 1;
	
		ncells++;

	    }

	}

    }



    

    // 3D

    else {


	cells.resize( (nx-1)*(ny-1)*(nz-1) );

	for(uint i = 0 ; i < (nx-1)*(ny-1)*(nz-1) ; i++)
	    cells[i].resize(8,0);
	

	uint ncells(0);

	for( uint k = 0 ; k < (nz-1) ; k++ ) {
    
	    for( uint j = 0 ; j < (ny-1) ; j++ ) {
	
		for( uint i = 0 ; i < (nx-1) ; i++ ) {
	

		    cells[ncells][0] = i   +   j*nx             +   k*nx*ny;

		    cells[ncells][1] = i   +   j*nx + 1         +   k*nx*ny;

		    cells[ncells][2] = i   +   j*nx + nx        +   k*nx*ny;

		    cells[ncells][3] = i   +   j*nx + nx + 1    +   k*nx*ny;

		    cells[ncells][4] = i   +   j*nx             +   (k+1)*nx*ny;

		    cells[ncells][5] = i   +   j*nx + 1         +   (k+1)*nx*ny;

		    cells[ncells][6] = i   +   j*nx + nx        +   (k+1)*nx*ny;

		    cells[ncells][7] = i   +   j*nx + nx + 1    +   (k+1)*nx*ny;
	
		    ncells++;


		}

	    }

	}
	


	
    }


}
