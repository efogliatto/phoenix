#ifndef PSEUDOPOTFORCE_H
#define PSEUDOPOTFORCE_H

#include <bForces.H>

#include <intForce.H>

#include <externalForce.H>


class pseudoPotForce {


protected:


    /** Interaction force */

    interactionForce* _Fi;


    /** Buoyant force */

    buoyantForce* _Fb;


    /** External force */

    externalForce _Fe;


    /** Reference to external density */

    const scalarField& _rho;


public:


    /** Constructor */

    pseudoPotForce( const std::string& dictName, const std::string& eqName, const latticeMesh& mesh, timeOptions& Time, const scalarField& rho );


    /** Destructor */

    ~pseudoPotForce();


    /** Total force  at node*/

    const std::vector<scalar> total( const uint& id ) const;


    /** Interaction force at node */

    const std::vector<scalar>& interaction( const uint& id ) const;

    const scalar interaction( const uint& i, const uint& j ) const { return _Fi->force(i,j); }

    void total( scalar Ft[3], const uint& id );


    /** Interaction force at node */

    const std::vector<scalar> buoyant( const uint& id ) const;


    /** Update forces */

    void update( scalarField& rho, scalarField& T );


    /** Sync force fields */

    void sync();


    /** Interaction potential */

    const scalar potential(const scalar& rho, const scalar& T, const scalar& cs2) const;
    
    
};


#endif // PSEUDOPOTFORCE_H
