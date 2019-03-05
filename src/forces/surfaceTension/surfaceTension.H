#ifndef SURFACETENSION_H
#define SURFACETENSION_H

#include <interactionForce.H>

#include <latticeMesh.H>

#include <scalarField.H>


class surfaceTension {

protected:


    /** Reference to lattice mesh */

    const latticeMesh& _mesh;


public:

    
    /** Surface tension model type */

    enum class stType { none, liST };
    
    
    /** Constructor */

    surfaceTension( const std::string& dictName, const std::string& eqName, const latticeMesh& mesh );


    /** Destructor */

    ~surfaceTension();    


    /** Surface Tension model */

    virtual const stType sTModel() const = 0;
    

    /** Surface tension at specific node */

    virtual const void ST( const uint& i, const scalarField& rho, const scalarField& T, std::vector<scalar>& C, interactionForce* _fi, const std::vector<scalar>& Tau ) const = 0;

    

};


#endif // SURFACETENSION_H
