#include <periodicBoundaryCorrection.H>

#include <iostream>

using namespace std;

void periodicBoundaryCorrection( vector< vector<int> >& nb,
				 const vector<periodicBnds>& periodicPairs,
				 const vector< vector<uint> >& points,
				 const unordered_map< string, vector<uint> >& boundaries ) {


    // Move over pairs and check

    for( const auto& bdpair : periodicPairs ) {


	
	// Indices list for each boundary

	vector<uint> bd1Points = boundaries.at( bdpair.bd1 );

	vector<uint> bd2Points = boundaries.at( bdpair.bd2 );


	
	// Move over point lists and check separation

	for( const auto& pt1 : bd1Points ) {

	    for( const auto& pt2 : bd2Points ) {


		// Separation vector between points
		
		int sep[3];

		for( uint j = 0 ; j < 3 ; j++ )
		    sep[j] = points[pt2][j] - points[pt1][j];
		


		if(  ( sep[0] == bdpair.direction[0] )   &&   ( sep[1] == bdpair.direction[1] )   &&   ( sep[2] == bdpair.direction[2] )  ) {

		    for( uint k = 0 ; k < nb[pt1].size() ; k++ ) {
			
			if( nb[pt1][k] == -1 )		    
			    nb[pt1][k] = nb[pt2][k];

			if( nb[pt2][k] == -1 )
			    nb[pt2][k] = nb[pt1][k];			

		    }

		}

		

	    }

	}
	

    }
    
    

}
