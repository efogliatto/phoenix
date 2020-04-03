#ifndef SINGLERANGEWITHCONTACT_H
#define SINGLERANGEWITHCONTACT_H

#include <singleRangeIntForce.H>


class singleRangeWithContact : public singleRangeIntForce {


protected:

    
    // Contact angle schemes

    enum class contactType{ exp, imp };

    contactType _scheme;


    // Limit angle

    scalar _limitAngle;
    

public:


    /** Constructor */

    singleRangeWithContact( const std::string& dictName, const std::string& eqName, const latticeMesh& mesh, timeOptions& Time );


    /** Destructor */

    ~singleRangeWithContact();


    /** Update force field */

    virtual void update( scalarField& rho, scalarField& T );    
    

};

#endif // SINGLERANGEWITHCONTACT_H
