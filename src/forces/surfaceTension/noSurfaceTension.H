#ifndef NOSURFACETENSION_H
#define NOSURFACETENSION_H

#include <surfaceTension.H>

class noSurfaceTension : public surfaceTension {

    
public:

    
    /** Constructor */

    noSurfaceTension( const std::string& dictName, const std::string& eqName, const latticeMesh& mesh );


    /** Destructor */

    ~noSurfaceTension();


    /** Surface Tension model */

    const stType sTModel() const { return stType::none; };
    

    /** Surface tension term at specific node */

    void ST( const uint& i, const scalarField& rho, const scalarField& T, std::vector<scalar>& C, interactionForce* _fi, const std::vector<scalar>& Tau );

    
};

#endif // NOSURFACETENSION_H
