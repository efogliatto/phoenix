#ifndef LISURFACETENSION_H
#define LISURFACETENSION_H

#include <surfaceTension.H>

class liSurfaceTension : public surfaceTension {


private:

    /** Model constant */

    scalar _kappa;

    
public:

    /** Constructor */

    liSurfaceTension( const std::string& dictName, const std::string& eqName, const latticeMesh& mesh );


    /** Destructor */

    ~liSurfaceTension();


    /** Additional surface term at specific node */

    const void ST( const uint& i, const scalarField& rho, const scalarField& T, std::vector<scalar>& C, interactionForce* _fi, const std::vector<scalar>& Tau ) const;
    

};

#endif // LISURFACETENSION_H
