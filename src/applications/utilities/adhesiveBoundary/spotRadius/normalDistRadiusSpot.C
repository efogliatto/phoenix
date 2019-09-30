#include <normalDistRadiusSpot.H>

#include <random>


using namespace std;


/** Constructor */
    
normalDistRadiusSpot::normalDistRadiusSpot() {}


/** Destructor */

normalDistRadiusSpot::~normalDistRadiusSpot() {}


/** Compute radius */

vector<uint> normalDistRadiusSpot::radius( const uint nspots, const uint mean, const uint dev ) const {

    vector<uint> rad(nspots);

    std::default_random_engine generator;
    
    std::normal_distribution<scalar> distribution( (scalar)mean,(scalar)dev );

    for( uint i = 0 ; i < nspots ; i++ )
	rad[i] = (uint) distribution(generator);

    
    return rad;

}
