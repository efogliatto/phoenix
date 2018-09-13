#ifndef INTFORCE_H
#define INTFORCE_H

#include <singleRangeIntForce.H>


class intForce {

public:

    /** Constructor */
    
    intForce( const std::string& dictName, const std::string& eqName, const latticeMesh& mesh, timeOptions& Time );


    /** Destructor */

    ~intForce();


    /** Overloaded () operator */

    interactionForce* operator() (const std::string& dictName, const std::string& eqName, const latticeMesh& mesh, timeOptions& Time);

};

#endif // INTFORCE_H
