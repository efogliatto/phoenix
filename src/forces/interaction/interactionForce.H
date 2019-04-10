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


    /** Flag: compute on boundary nodes */

    bool _computeOnBnd;


    /** Interaction weights */

    std::vector<scalar> _weights;


    /** Geometric contact angle */

    std::map< uint, scalar > _contactAngle;

    std::map< uint, std::vector<scalar> > _hysteresis;

    
    /** Flag for geometric contact line location*/

    bool _withGeomContact;



    /** Read geometric contact angle properties */

    const void readGeometricContact( const std::string& eqname );

    const scalar apparentContactAngle( const scalarField& rho, const std::string& bname );    
    

    

public:


    /** Constructor */

    interactionForce( const std::string& dictName, const std::string& eqName, const latticeMesh& mesh, timeOptions& Time );


    /** Destructor */

    ~interactionForce();


    /** Force at specific node */

    const std::vector<scalar>& force( const uint& i ) { return _force[i]; }

    const scalar force( const uint& i, const uint& j ) { return _force[i][j]; }


    /** Set force at specific node */

    const void set( const uint& i, const std::vector<scalar>& Fint );


    /** Update force field */

    virtual void update( scalarField& rho, scalarField& T ) = 0;


    /** Interaction potential */

    const scalar potential( const scalar& rho, const scalar& T, const scalar& cs2 ) const;


    /** Adequately signed potential strength */

    const scalar signedPotentialStrength( const scalar& rho, const scalar& T, const scalar& cs2 ) const;


    /** Synv force field */

    const void sync() { _force.sync(); }


    /** Potential weights */

    const std::vector<scalar>& weights() const { return _weights; }
    

};

#endif // INTERACTIONFORCE_H
