#ifndef NORMALDISTRADIUSSPOT_H
#define NORMALDISTRADIUSSPOT_H

#include <spotRadius.H>

class normalDistRadiusSpot : public spotRadius {

public:


    /** Constructor */
    
    normalDistRadiusSpot();


    /** Destructor */

    ~normalDistRadiusSpot();


    /** Compute radius */

    std::vector<uint> radius( const uint nspots, const uint mean, const uint dev ) const;

};

#endif // NORMALDISTRADIUSSPOT_H
