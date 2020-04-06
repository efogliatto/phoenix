#ifndef SINGLERANGEMIXEDWITHCONTACT_H
#define SINGLERANGEMIXEDWITHCONTACT_H

#include <singleRangeMixedIntForce.H>


class singleRangeMixedWithContact : public singleRangeMixedIntForce {


protected:

    
    // Contact angle schemes

    enum class contactType{ exp, imp };

    contactType _scheme;


    // Limit angle

    scalar _limitAngle;
    

public:


    /** Constructor */

    singleRangeMixedWithContact( const std::string& dictName, const std::string& eqName, const latticeMesh& mesh, timeOptions& Time );


    /** Destructor */

    ~singleRangeMixedWithContact();


    /** Update force field */

    virtual void update( scalarField& rho, scalarField& T );    
    

};

#endif // SINGLERANGEMIXEDWITHCONTACT_H
