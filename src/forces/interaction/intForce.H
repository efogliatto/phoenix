#ifndef INTFORCE_H
#define INTFORCE_H

#include <singleRangeMixedIntForce.H>


class intForce {

public:

    
    /** Force creator */
    
    interactionForce* create( const std::string& dictName, const std::string& eqName, const latticeMesh& mesh, timeOptions& Time );


};

#endif // INTFORCE_H
