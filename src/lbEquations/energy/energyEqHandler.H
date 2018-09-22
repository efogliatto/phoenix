#ifndef ENERGYEQHANDLER
#define ENERGYEQHANDLER

#include <EEquation.H>

#include <energyBndCreator.H>


class energyEqHandler {

protected:

    
    /** LB energy equation */

    energyEquation* _equation;


    /** Boundary conditions */

    std::vector< energyBndCond* > _boundaries;



public:

    /** Constructor */

    energyEqHandler(const std::string& name,
		    const latticeMesh& mesh_,
		    timeOptions& Time_,
		    pdfField& pdf_,
		    const scalarField& rho_,
		    const vectorField& U_,
		    scalarField& T_);


    /** Destructor */

    ~energyEqHandler();


    /** Collision */

    const void collision() { _equation->collision(); }


    /** Streaming */

    const void streaming() { _equation->streaming(); }


    /** Update macroscopic temperature */

    const void updateMacroTemperature() { _equation->updateMacroTemperature(); }
    

};

#endif // ENERGYEQHANDLER
