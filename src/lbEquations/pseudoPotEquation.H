#ifndef PSEUDOPOTEQUATION_H
#define PSEUDOPOTEQUATION_H

#include <lbEquation.H>

#include <EOSCreator.H>


class pseudoPotEquation : public lbEquation {


protected:

    
    /** Reference to external density field */

    scalarField& rho;


    /** Reference to external velocity field */

    vectorField& U;


    /** Reference to external temperature field */

    const scalarField& T;


    /** Interaction force */

    vectorField Fi;


    /** Total force */

    vectorField Ft;


    /** Equation of state */

    EOS* eos;
    


    

    /** Compute local density */

    scalar localDensity( const uint& id );


    /** Compute local velocity */
    /** Density MUST be already updated */

    const void localVelocity( std::vector<scalar>& v, const uint& id );    


    



public:

    
    /** Default constructor */

    pseudoPotEquation( const std::string& name,
		       const latticeMesh& mesh_,
		       timeOptions& Time_,
		       pdfField pdf_,
		       scalarField& rho_,
		       vectorField& U_,
		       scalarField& T_);


    /** Default destructor */

    ~pseudoPotEquation();



    /** Collision process */

    virtual const void collision() = 0;


    /** Update macroscopic density */

    const void updateMacroDensity();

    
    /** Update macroscopic velocity */

    const void updateMacroVelocity();


    /** Set pdf to equilibrium values */

    virtual const void setEquilibrium() = 0;


    
};


#endif // PSEUDOPOTEQUATION_H
