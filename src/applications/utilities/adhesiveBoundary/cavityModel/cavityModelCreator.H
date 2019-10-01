#ifndef CAVITYMODELCREATOR_H
#define CAVITYMODELCREATOR_H

#include <flatCavity.H>

#include <linearCavity.H>

#include <map>

#include <memory>


class cavityModelCreator {

    
protected:

    /** Model types */

    enum class cavityType{flat, linear};

    std::map<std::string, cavityType> _cavityMapType;


public:

    /** Model creation */
    
    std::unique_ptr<cavityModel> create( const std::string& bdname );   

    
};


#endif // CAVITYMODELCREATOR_H
