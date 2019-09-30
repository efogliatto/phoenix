#ifndef SPOTSAMPLECREATOR_H
#define SPOTSAMPLECREATOR_H

#include <uniformSpots.H>

#include <randomSpots.H>

#include <map>

#include <memory>


class spotSampleCreator {

    
protected:

    /** Model types */

    enum class spotType{uniform, random};

    std::map<std::string, spotType> _spotMapType;


public:

    /** Model creation */
    
    std::unique_ptr<spotSample> create( const std::string& sptype, const std::string& bdname );   

    
};

#endif // SPOTSAMPLECREATOR_H
