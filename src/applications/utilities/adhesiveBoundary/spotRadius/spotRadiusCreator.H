#ifndef SPOTRADIUSCREATOR_H
#define SPOTRADIUSCREATOR_H

#include <fixedRadiusSpot.H>

#include <normalDistRadiusSpot.H>

#include <map>

#include <memory>


class spotRadiusCreator {

    
protected:

    
    /** Model types */

    enum class radiusType{fixed, normal};

    std::map<std::string, radiusType> _radiusMapType;


public:

    /** Model creation */
    
    std::unique_ptr<spotRadius> create( const std::string& rtype );   

    
};

#endif // SPOTRADIUSCREATOR_H
