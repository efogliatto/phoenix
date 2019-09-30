#include <fixedRadiusSpot.H>

#include <algorithm>

using namespace std;


/** Constructor */
    
fixedRadiusSpot::fixedRadiusSpot() {}


/** Destructor */

fixedRadiusSpot::~fixedRadiusSpot() {}


/** Compute radius */

vector<uint> fixedRadiusSpot::radius( const uint nspots, const uint mean, const uint dev ) const {

    vector<uint> rad(nspots);

    std::fill( rad.begin(), rad.end(), mean );
    
    return rad;

}
