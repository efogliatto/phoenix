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




    // Index mapping for new points

    vector<int> oldToNew;

    oldToNew.resize( basePoints.size(), -1 );


    uint nPoints(0);

    for( uint i = 0 ; i < baseCells.size() ; i++ ) {

	if( isInside[i] ) {

	    for( uint j = 0 ; j < baseCells[i].size() ; j++ ) {

		uint pt = baseCells[i][j];

		if( oldToNew[pt] == -1 ) {

		    oldToNew[pt] = nPoints;
		    
		    nPoints++;

		}
		
	    }

	}

    }



    // Resize new point array and copy

    newPoints.resize( nPoints );

    for( uint i = 0 ; i < nPoints ; i++ )
	newPoints[i].resize(3,0);


    for( uint i = 0 ; i < basePoints.size() ; i++ ) {

	int pt = oldToNew[i];

	if( pt != -1 ) {

	    for(uint j = 0 ; j < 3 ; j++)
		newPoints[pt][j] = basePoints[i][j];

	}
	    

    }



    // Resize new cell array and copy

    newCells.resize( nCells );

    for( uint i = 0 ; i < nCells ; i++ )
	newCells[i].resize( baseCells[0].size() ,0);

    {

	uint k(0);

	for( uint i = 0 ; i < baseCells.size() ; i++ ) {

	    if( isInside[i] ) {

		for( uint j = 0 ; j < baseCells[i].size() ; j++ )
		    newCells[k][j] = oldToNew[ baseCells[i][j] ];

		k++;

	    }

	}

    }
    

    

}
