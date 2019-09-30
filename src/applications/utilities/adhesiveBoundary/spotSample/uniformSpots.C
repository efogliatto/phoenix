#include <uniformSpots.H>


using namespace std;


/** Constructor */

uniformSpots::uniformSpots(const std::string& bdname) : spotSample(bdname) {

    
    // Distribute locations uniformly according node indices

    uint step = (uint)( _nodes.size() / ( _spots.size()+1 ) );

    for(uint i = 0 ; i < _spots.size() ; i++)
    	_spots[i].first = _nodes[ (i+1)*step ];

}


/** Destructor */

uniformSpots::~uniformSpots(){}
