#ifndef SPOTRADIUS_H
#define SPOTRADIUS_H

#include <dataTypes.H>

#include <dictionary.H>


class spotRadius {

public:


    /** Constructor */
    
    spotRadius();


    /** Destructor */

    ~spotRadius();


    /** Compute radius */

    virtual std::vector<uint> radius( const uint nspots, const uint mean, const uint dev ) const = 0;

};

#endif // SPOTRADIUS_H
