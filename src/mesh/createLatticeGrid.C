#include <createLatticeGrid.H>

using namespace std;

void createLatticeGrid( vector< vector<uint> >& points, const latticeModel* lbmodel, const uint nx, const uint ny, const uint nz) {


    
    // Create grid according to lattice model


    // 2D
    
    if( lbmodel->d() == 2 ) {

	points.resize( nx * ny );

	for( uint i = 0 ; i < nx*ny ; i++ )
	    points[i].resize(3,0);


	for( uint j = 0 ; j < ny ; j++) {
	
	    for( uint i = 0 ; i < nx ; i++) {
			    
		points[i+j*nx][0] = i;
		points[i+j*nx][1] = j;		
		
	    }
	}
	
	

    }


    // 3D

    else {


	points.resize( nx * ny * nz );

	for( uint i = 0 ; i < nx*ny*nz ; i++ )
	    points[i].resize(3,0);


	for( uint k = 0 ; k < nz ; k++) {
    
	    for( uint j = 0 ; j < ny ; j++) {
	
		for( uint i = 0 ; i < nx ; i++) {

		    uint idx = i + j*nx + k*nx*ny;
			    
		    points[idx][0] = i;
		
		    points[idx][1] = j;

		    points[idx][2] = k;
		
		}
	    }

	}	
	

    }


}
