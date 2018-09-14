#ifndef PSEUDOPOTFORCE_H
#define PSEUDOPOTFORCE_H

#include <bForces.H>

#include <intForce.H>


class pseudoPotForce {


protected:


    /** Interaction force */

    interactionForce* _Fi;


    /** Buoyant force */

    buoyantForce* _Fb;


public:


    /** Constructor */

    pseudoPotForce( const std::string& dictName, const std::string& eqName, const latticeMesh& mesh, timeOptions& Time );


    /** Destructor */

    ~pseudoPotForce();


    /** Total force  at node*/

    const std::vector<scalar> total( const uint& id ) const;


    /** Interaction force at node */

    const std::vector<scalar> interaction( const uint& id ) const;


    /** Interaction force at node */

    const std::vector<scalar> buoyant( const uint& id ) const;


    /** Update forces */

    void update( scalarField& rho, scalarField& T );


    /** Sync force fields */

    void sync();
    
    
};


#endif // PSEUDOPOTFORCE_H
