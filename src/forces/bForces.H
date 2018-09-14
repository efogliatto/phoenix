#ifndef BFORCES_H
#define BFORCES_H

#include <fixedDensityBForce.H>

class bForces {

public:

    buoyantForce* create( const std::string& dictName, const std::string& eqName );

};

#endif // BFORCES_H
