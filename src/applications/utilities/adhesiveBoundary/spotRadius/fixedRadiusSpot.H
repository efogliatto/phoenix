#ifndef FIXEDRADIUSSPOT_H
#define FIXEDRADIUSSPOT_H

#include <spotRadius.H>

class fixedRadiusSpot : public spotRadius {

public:


    /** Constructor */
    
    fixedRadiusSpot();


    /** Destructor */

    ~fixedRadiusSpot();


    /** Compute radius */

    std::vector<uint> radius( const uint nspots, const uint mean, const uint dev ) const;

};

#endif // FIXEDRADIUSSPOT_H
