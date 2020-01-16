#ifndef SPOTSAMPLE_H
#define SPOTSAMPLE_H

#include <dataTypes.H>

#include <string>

#include <vector>

#include <cavityModelCreator.H>


class spotSample {
    
    
protected:

    
    /** Boundary indices */

    std::vector<uint> _nodes;


    /** Nodes location */

    std::vector< std::vector<uint> > _location;
    

    /** Spots. Centre (node id) and radius */

    std::vector< std::pair<uint,uint> > _spots;


    /** Cavity model */

    std::unique_ptr<cavityModel> _cavity;



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


    /** Compute adhesive coefficients */

    virtual const std::vector< std::pair<uint, scalar> > computeAds() const;

    

};

#endif // SPOTSAMPLE_H
