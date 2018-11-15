#ifndef ADSORPTIONFORCE_H
#define ADSORPTIONFORCE_H

#include <EOSCreator.H>

#include <latticeMesh.H>


class adsorptionForce {

protected:


    /** Reference to lattice mesh */

    const latticeMesh& _mesh;


    /** Equation of state */

    EOS* eos;


public:

    
    /** Constructor */

    adsorptionForce( const std::string& dictName, const std::string& eqName, const latticeMesh& mesh );


    /** Destructor */

    ~adsorptionForce();


    /** Force at specific node */

    virtual const void force( const uint& i, scalar f[3] ) const = 0;

    

};


#endif // ADSORPTIONFORCE_H
