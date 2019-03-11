#ifndef PSEUDOPOTFORCE_H
#define PSEUDOPOTFORCE_H

#include <bForces.H>

#include <intForce.H>

#include <adsForceCreator.H>

#include <externalForce.H>

#include <stCreator.H>


class pseudoPotForce {


protected:


    /** Interaction force */

    interactionForce* _Fi;


    /** Buoyant force */

    buoyantForce* _Fb;


    /** Adhesive force */

    adhesiveForce* _Fads;


    /** External force */

    externalForce _Fe;


    /** Additional surface tension term */

    surfaceTension* _St;

    

    /** Reference to external density */

    const scalarField& _rho;


    /** Reference to external temperature */

    const scalarField& _T;


    /** Reference to latticeMesh */

    const latticeMesh& _mesh;
    


public:


    /** Constructor */

    pseudoPotForce( const std::string& dictName,
		    const std::string& eqName,
		    const latticeMesh& mesh,
		    timeOptions& Time,
		    const scalarField& rho,
		    const scalarField& T );


    /** Destructor */

    ~pseudoPotForce();


    /** Total force  at node*/

    const std::vector<scalar> total( const uint& id ) const;


    /** Interaction force at node */

    const std::vector<scalar>& interaction( const uint& id ) const;

    const scalar interaction( const uint& i, const uint& j ) const { return _Fi->force(i,j); }

    void total( scalar Ft[3], const uint& id ) const;


    /** Set interaction force at node */

    const void setInteraction( const uint& i, const std::vector<scalar>& Fint ) { _Fi->set(i, Fint); }


    /** Interaction force at node */

    const std::vector<scalar> buoyant( const uint& id ) const;


    /** Update forces */

    void update( scalarField& rho, scalarField& T );


    /** Sync force fields */

    void sync();


    /** Interaction potential */

    const scalar potential(const scalar& rho, const scalar& T, const scalar& cs2) const;


    /** Adequately signed potential strength */

    const scalar signedPotential( const scalar& rho, const scalar& T, const scalar& cs2 ) const;


    /** Additional surface tension term at node */

    const void addSurfaceTension( const uint& i, std::vector<scalar>& C, const std::vector<scalar>& Tau ) const;


    /** Surface tension model */

    const surfaceTension::stType sTModel() const { return _St->sTModel(); }
    
    
};


#endif // PSEUDOPOTFORCE_H
