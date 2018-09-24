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


    // /** Boundary conditions */

    // std::vector< ppBndCond* > _boundaries;


  
    



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


    // /** Update boundaries */

    // void updateBoundaries();




    /** Compute local density */

    const scalar localDensity( const uint& id ) const;


    /** Compute local velocity */

    const void localVelocity( std::vector<scalar>& v, const uint& id, const bool updDens = false ) const;

    
    
};


#endif // PSEUDOPOTEQUATION_H
