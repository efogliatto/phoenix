#ifndef PHIBASEDMOD_H
#define PHIBASEDMOD_H

#include <adhesiveForce.H>

class phiBasedMod : public adhesiveForce {

    
public:

    
    /** Constructor */

    phiBasedMod( const std::string& dictName, const std::string& eqName, const latticeMesh& mesh, const interactionForce* Fi, timeOptions& Time );


    /** Destructor */

    ~phiBasedMod();



    /** Update force field */

    void update( scalarField& rho, scalarField& T );

    
};

#endif // PHIBASEDMOD_H
