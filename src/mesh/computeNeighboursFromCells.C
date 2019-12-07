#include <computeNeighboursFromCells.H>

using namespace std;

void computeNeighboursFromCells( vector< vector<int> >& nb, vector< vector<uint> >& points, vector< vector<uint> >& cells, const latticeModel* lbmodel ) {


    // Lattice constants

    const uint nPoints( points.size() );

    const uint nCells( cells.size() );

    const uint q = lbmodel->q();

    const uint ppc( cells[0].size() );

    const vector<uint> reverse = lbmodel->reverse();

    
    
    // Resize neighbour array

    nb.resize( nPoints );

    for( uint i = 0 ; i < nPoints ; i++ )
	nb[i].resize( q, -1 );


    
    // Move over cells and check neighbouring

    int sep[3];

    for( uint i = 0 ; i < nCells ; i++ ) {


	// Check for every point in cell

	for( uint j = 0 ; j < ppc ; j++ ) {

	    for( uint k = j ; k < ppc ; k++ ) {


		// Compute separation between points k and j, and check if it matches any lattice velocity
		
		sep[0] = points[ cells[i][k] ][0] - points[ cells[i][j] ][0];
		sep[1] = points[ cells[i][k] ][1] - points[ cells[i][j] ][1];
		sep[2] = points[ cells[i][k] ][2] - points[ cells[i][j] ][2];


		int vid = lbmodel->velocityIndex( sep[0], sep[1], sep[2] );

		if(vid != -1) {

		    nb[ cells[i][j] ][ reverse[vid] ] = cells[i][k];

		    nb[ cells[i][k] ][ vid ] = cells[i][j];		    

		}
		    

	    }	    

	}
	

    }

}
