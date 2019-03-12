#ifndef ADHESIVEFORCE_H
#define ADHESIVEFORCE_H

#include <EOSCreator.H>

#include <latticeMesh.H>

#include <scalarField.H>

#include <vectorField.H>

#include <intForce.H>


/*

  Adhesive force applies only on local nodes (like interaction force)

 */


class adhesiveForce {

protected:


    /** Reference to lattice mesh */

    const latticeMesh& _mesh;


    /** Reference to related interaction force */

    const interactionForce* _Fi;


    /** Force field */

    vectorField _force;

    
    /** Map of nodes where adhesive force applies */

    std::map<uint, scalar> _Gads;
    


public:

    
    /** Constructor */

    adhesiveForce( const std::string& dictName, const std::string& eqName, const latticeMesh& mesh, const interactionForce* Fi, timeOptions& Time );


    /** Destructor */

    ~adhesiveForce();

    
    /** Update force field */

    virtual void update( scalarField& rho, scalarField& T ) = 0;


    /** Force at specific node */

    const std::vector<scalar>& force( const uint& i ) { return _force[i]; }

    const scalar force( const uint& i, const uint& j ) { return _force[i][j]; }


    /** Sync force field */

    const void sync() { _force.sync(); }
    

    

};


#endif // ADHESIVEFORCE_H
