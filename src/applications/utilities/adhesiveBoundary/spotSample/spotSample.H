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
    

    /** Spots */

    std::vector< std::pair<uint,uint> > _spots;



    /** Read nodes on boundary */

    const void readNodes( const std::string& bdname );


    /** Read position of nodes on boundary */

    const void readLocations( const std::string& bdname );


    /** Create spots */

    const void createSpots( const std::string& bdname );

    

public:


    /** Default constructor */

    spotSample(const std::string& bdname);


    /** Destructor */

    ~spotSample();

    

};

#endif // SPOTSAMPLE_H
