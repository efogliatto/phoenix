#ifndef PSEUDOPOTEQHANDLER
#define PSEUDOPOTEQHANDLER

#include <PPEquation.H>

#include <ppBndCreator.H>


class pseudoPotEqHandler {

protected:

    
    /** LB pseudopotential equation */

    pseudoPotEquation* _equation;


    /** Boundary conditions */

    std::vector< ppBndCond* > _boundaries;



public:

    /** Constructor */

    pseudoPotEqHandler( const std::string& name,
			const latticeMesh& mesh_,
			timeOptions& Time_,
			pdfField& pdf_,
			scalarField& rho_,
			vectorField& U_,
			scalarField& T_ );


    /** Destructor */

    ~pseudoPotEqHandler();


    /** Collision */

    const void collision() { _equation->collision(); }


    /** Streaming */

    const void streaming() { _equation->streaming(); }


    /** Update macroscopic density */

    const void updateMacroDensity() { _equation->updateMacroDensity(); }

    
    /** Update macroscopic velocity */

    const void updateMacroVelocity() { _equation->updateMacroVelocity(); }
    

    /** Update boundaries */

    const void updateBoundaries();
    

};

#endif // PSEUDOPOTEQHANDLER
