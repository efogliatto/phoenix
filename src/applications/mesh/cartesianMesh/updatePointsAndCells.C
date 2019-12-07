#include <updatePointsAndCells.H>

using namespace std;


void updatePointsAndCells( vector< vector<uint> >& basePoints,
			   vector< vector<uint> >& baseCells,
			   vector< vector<uint> >& newPoints,
			   vector< vector<uint> >& newCells,
			   const vector<bool>& isInside ) {


    // Count valid cells

    uint nCells(0);
    
    for( auto in : isInside ) {

	if( in )
	    nCells++;

    }


    newCells.resize( nCells );



    // Index mapping for new points

    vector<int> oldToNew;

    oldToNew.resize( basePoints.size(), -1 );


    uint nPoints(0);

    for( uint i = 0 ; i < baseCells.size() ; i++ ) {

	if( isInside[i] ) {

	    for( uint j = 0 ; j < baseCells[i].size() ; j++ ) {

		

	    }

	}

    }

    

}
