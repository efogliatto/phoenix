#include <randomSpots.H>

#include <random>


using namespace std;


/** Constructor */

randomSpots::randomSpots(const std::string& bdname) : spotSample(bdname) {

    
    // Distribute locations randomly according node indices

    std::random_device rd;
    
    std::mt19937 eng(rd()); 

    std::uniform_int_distribution<> distr(0, _nodes.size()-1); 
    

    for(uint i = 0 ; i < _spots.size() ; i++)
    	_spots[i].first = distr(eng);

}


/** Destructor */

randomSpots::~randomSpots(){}
