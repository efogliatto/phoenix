#ifndef RNDSPOTS_H
#define RNDSPOTS_H

#include <phiBasedMod.H>

class rndSpots : public phiBasedMod {

    
public:

    
    /** Constructor */

    rndSpots( const std::string& dictName, const std::string& eqName, const latticeMesh& mesh, const interactionForce* Fi, timeOptions& Time );


    /** Destructor */

    ~rndSpots();



    /** Update force field */

    void update( scalarField& rho, scalarField& T );

    
};

#endif // RNDSPOTS_H
