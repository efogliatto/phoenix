#ifndef INTFORCE_H
#define INTFORCE_H

#include <singleRangeMixedWithContact.H>

#include <singleRangeWithContact.H>


class intForce {

    
protected:

    
    // Interaction types

    enum class ifType{ singleRange, singleRangeMixed, singleRangeContact, singleRangeMixedContact };

    std::map<std::string, ifType> _ifMapType;

    
public:

    
    /** Force creator */
    
    interactionForce* create( const std::string& dictName, const std::string& eqName, const latticeMesh& mesh, timeOptions& Time );


};

#endif // INTFORCE_H
