#ifndef NOSURFACETENSION_H
#define NOSURFACETENSION_H

#include <surfaceTension.H>

class noSurfaceTension : public surfaceTension {

    
public:

    
    /** Constructor */

    noSurfaceTension( const std::string& dictName, const std::string& eqName, const latticeMesh& mesh );


    /** Destructor */

    ~noSurfaceTension();


    /** Surface tension term at specific node */

    const void ST( const uint& i, const scalarField& rho, const scalarField& T, std::vector<scalar>& C ) const;

    
};

#endif // NOSURFACETENSION_H
