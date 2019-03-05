#ifndef SURFACETENSION_H
#define SURFACETENSION_H

// #include <EOSCreator.H>

#include <latticeMesh.H>

#include <scalarField.H>


class surfaceTension {

protected:


    /** Reference to lattice mesh */

    const latticeMesh& _mesh;


    // /** Equation of state */

    // EOS* eos;


public:

    
    /** Constructor */

    surfaceTension( const std::string& dictName, const std::string& eqName, const latticeMesh& mesh );


    /** Destructor */

    ~surfaceTension();


    /** Surface tension at specific node */

    virtual const void ST( const uint& i, const scalarField& rho, const scalarField& T, std::vector<scalar>& C ) const = 0;

    

};


#endif // SURFACETENSION_H
