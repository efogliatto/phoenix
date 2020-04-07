#ifndef LISURFACETENSIONCONTACT_H
#define LISURFACETENSIONCONTACT_H

#include <liSurfaceTension.H>

class liSurfaceTensionContact : public liSurfaceTension {


private:

    /** Geometric contact angle */

    std::map< uint, scalar > _contactAngle;


    /** Hysteresis range */
    
    std::map< uint, std::vector<scalar> > _hysteresis;


    /** Limit angle */

    scalar _limitAngle;
    
    
    /** Read geometric contact angle properties */

    const void readGeometricContact( const std::string& eqname );

    
public:

    /** Constructor */

    liSurfaceTensionContact( const std::string& dictName, const std::string& eqName, const latticeMesh& mesh );


    /** Destructor */

    ~liSurfaceTensionContact();

    
    /** Additional surface term at specific node */

    void ST( const uint& i, const scalarField& rho, const scalarField& T, std::vector<scalar>& C, interactionForce* _fi, const std::vector<scalar>& Tau );   

};

#endif // LISURFACETENSIONCONTACT_H
