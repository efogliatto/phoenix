#ifndef SPOTSAMPLE_H
#define SPOTSAMPLE_H

#include <dataTypes.H>

#include <string>

#include <vector>


class spotSample {

protected:

    
    /** Boundary indices */

    std::vector<uint> _nodes;


    /** Nodes location */

    std::vector< std::vector<uint> > _location;
    

    /** Total number of spots */
    
    uint _nSpots;


    /** Mean radius */

    scalar _radius;
    
    

public:


    /** Default constructor */

    spotSample(const std::string& bdname);


    /** Destructor */

    ~spotSample();

    

};

#endif // SPOTSAMPLE_H
