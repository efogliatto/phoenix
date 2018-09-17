#ifndef ENERGYEQUATION_H
#define ENERGYEQUATION_H

#include <lbEquation.H>



class energyEquation : public lbEquation {


protected:

    
    /** Reference to external density field */

    const scalarField& rho;


    /** Reference to external velocity field */

    const vectorField& U;


    /** Reference to external temperature field */

    scalarField& T;


    // /** Boundary conditions */

    // std::vector< ppBndCond* > _boundaries;


    

    /** Compute local temperature */

    scalar localTemperature( const uint& id );



public:

    
    /** Default constructor */

    energyEquation( const std::string& name,
		       const latticeMesh& mesh_,
		       timeOptions& Time_,
		       pdfField& pdf_,
		       const scalarField& rho_,
		       const vectorField& U_,
		       scalarField& T_);


    /** Default destructor */

    ~energyEquation();



    /** Collision process */

    virtual const void collision() = 0;


    /** Update macroscopic temperature */

    const void updateMacroTemperature();


    /** Set pdf to equilibrium values */

    virtual const void setEquilibrium() = 0;


    // /** Update boundaries */

    // void updateBoundaries();

    
};


#endif // ENERGYEQUATION_H
