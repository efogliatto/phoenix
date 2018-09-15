#ifndef INTERACTIONFORCE_H
#define INTERACTIONFORCE_H

#include <EOSCreator.H>

#include <vectorField.H>

#include <scalarField.H>



class interactionForce {


protected:


    /** Reference to lattice mesh */

    const latticeMesh& _mesh;
    
    
    /** Force field */

    vectorField _force;


    /** Equation of state */

    EOS* eos;


    /** Interaction strength */

    scalar _G;



public:


    /** Constructor */

    interactionForce( const std::string& dictName, const std::string& eqName, const latticeMesh& mesh, timeOptions& Time );


    /** Destructor */

    ~interactionForce();


    /** Force at specific node */

    const std::vector<scalar>& force( const uint& i ) { return _force[i]; }

    const scalar force( const uint& i, const uint& j ) { return _force[i][j]; }
    


    /** Update force field */

    virtual void update( scalarField& rho, scalarField& T ) = 0;


    /** Interaction potential */

    const scalar potential( const scalar& rho, const scalar& T, const scalar& cs2 ) const;


    /** Synv force field */

    const void sync() { _force.sync(); }
    

};

#endif // INTERACTIONFORCE_H
