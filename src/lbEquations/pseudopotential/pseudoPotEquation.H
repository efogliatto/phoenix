#ifndef PSEUDOPOTEQUATION_H
#define PSEUDOPOTEQUATION_H

#include <lbEquation.H>

#include <EOSCreator.H>

#include <pseudoPotForce.H>

// #include <ppBndCreator.H>



class pseudoPotEquation : public lbEquation {


protected:

    
    /** Reference to external density field */

    scalarField& rho;


    /** Reference to external velocity field */

    vectorField& U;


    /** Reference to external temperature field */

    scalarField& T;


    /** Model related forces */

    pseudoPotForce F;
 
    



public:

    
    /** Default constructor */

    pseudoPotEquation( const std::string& name,
		       const latticeMesh& mesh_,
		       timeOptions& Time_,
		       pdfField& pdf_,
		       scalarField& rho_,
		       vectorField& U_,
		       scalarField& T_);


    /** Default destructor */

    ~pseudoPotEquation();



    /** Collision process */

    virtual const void collision() = 0;


    /** Update macroscopic density */

    const void updateMacroDensity();

    const void updateMacroDensity(const uint& first, const uint& last);    

    
    /** Update macroscopic velocity */

    const void updateMacroVelocity();


    /** Set pdf to equilibrium values */

    virtual const void setEquilibrium() = 0;


    /** Forced interface for equilibrium distribution */

    virtual const void eqPS( std::vector<scalar>& n, const scalar& rho_, const std::vector<scalar>& U_ ) const;


    /** Access to interaction force */

    const scalar intForce(const uint& i, const uint& j) const { return F.interaction(i,j); }

    const void setIntForce(const uint& i, const std::vector<scalar>& Fint)  { F.setInteraction(i, Fint); }


    /** Access to total force at node */

    const void totalForce( scalar Ft[3], const uint& i ) const { F.total(Ft,i); }


    /** Update forces */

    const void updateForces() {  F.update(rho,T);  }



    /** Compute local density */

    const scalar localDensity( const uint& id ) const;


    
    /** Compute local density with unknowns. Useful in wet-node approach */

    const scalar localDensityWithUnknowns( const uint& id, latticeMesh::normalType& ntype ) const;    

    

    /** Compute local velocity */

    const void localVelocity( std::vector<scalar>& v, const uint& id, const bool updDens = false ) const;

    

    /** Update potential as scalar field */

    const void updatePotential( scalarField& phi );


    
    /** Compute and set pressure field */

    virtual const void pressure( const scalarField& phi, scalarField& p );
    
    
};


#endif // PSEUDOPOTEQUATION_H
