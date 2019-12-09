#include <periodicBoundaryCorrection.H>

using namespace std;

void periodicBoundaryCorrection( vector< vector<int> >& nb,
				 const map< pair<string,string>, vector<scalar> >& periodicPairs,
				 const vector< vector<uint> >& points,
				 const map< string, vector<uint> >& boundaries ) {


    // Move over pairs and check

    for( const auto& bdpair : periodicPairs ) {


	// Boundaries names

	string bd1_name = bdpair.first.first;

	string bd2_name = bdpair.first.second;

	vector<scalar> sepVec = bdpair.second;
	
	

	
	// Indices list for each boundary

	vector<uint> bd1 = boundaries.at( bd1_name );

	vector<uint> bd2 = boundaries.at( bd2_name );


	
	// Move over point lists and check separation

	for( const auto& pt1 : bd1 ) {

	    for( const auto& pt2 : bd2 ) {


		// Separation vector between points
		
		int sep[3];

		for( uint j = 0 ; j < 3 ; j++ )
		    sep[j] = points[pt2][j] - points[pt1][j];


		scalar sepMag(0);

		for( uint j = 0 ; j < 3 ; j++ )
		    sepMag += sep[j] * sep[j];

		sepMag = sepMag * sepMag;


		
		// Difference with forced direction

		scalar forcedSepMag(0);

		for( uint j = 0 ; j < 3 ; j++ )
		    forcedSepMag += sep[j] * sepVec[j];

		forcedSepMag = forcedSepMag * forcedSepMag;
		


		// Apply correction if magnitudes are equal
		
		if( sepMag == forcedSepMag ) {

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
