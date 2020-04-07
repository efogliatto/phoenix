#ifndef LISURFACETENSION_H
#define LISURFACETENSION_H

#include <surfaceTension.H>

class liSurfaceTension : public surfaceTension {


protected:

    /** Model constant */

    scalar _kappa;

    
public:

    /** Constructor */

    liSurfaceTension( const std::string& dictName, const std::string& eqName, const latticeMesh& mesh );


    /** Destructor */

    ~liSurfaceTension();


    /** Surface Tension model */

    virtual const stType sTModel() const { return stType::liST; };

    
    /** Additional surface term at specific node */

    virtual const void ST( const uint& i, const scalarField& rho, const scalarField& T, std::vector<scalar>& C, interactionForce* _fi, const std::vector<scalar>& Tau ) const;


    /** Model constant value */

    const scalar kappa() const;
    

};

#endif // LISURFACETENSION_H
