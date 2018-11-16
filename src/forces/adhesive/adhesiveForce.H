#ifndef ADHESIVEFORCE_H
#define ADHESIVEFORCE_H

#include <EOSCreator.H>

#include <latticeMesh.H>

#include <scalarField.H>


class adhesiveForce {

protected:


    /** Reference to lattice mesh */

    const latticeMesh& _mesh;


    /** Equation of state */

    EOS* eos;


public:

    
    /** Constructor */

    adhesiveForce( const std::string& dictName, const std::string& eqName, const latticeMesh& mesh );


    /** Destructor */

    ~adhesiveForce();


    /** Force at specific node */

    virtual const void force( const uint& i, const scalarField& rho, const scalarField& T, scalar f[3] ) const = 0;

    

};


#endif // ADHESIVEFORCE_H
