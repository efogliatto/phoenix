#ifndef PSEUDOPOTEQHANDLER
#define PSEUDOPOTEQHANDLER

#include <PPEquation.H>

#include <ppBndCreator.H>


class pseudoPotEqHandler {

protected:


    /** Mesh reference */

    const latticeMesh& _mesh;
    

    /** Density reference */

    const scalarField& _rho;
    
    
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

    const void updateMacroVelocity();
    

    /** Update boundaries */

    const void updateBoundaries();


    /** Update potential as scalar field  */

    const void updatePotential( scalarField& phi );


    /** Compute and set pressure field */

    const void pressure( const scalarField& phi, scalarField& p );


    // /** Check for contact points over surfaces */

    // std::vector<uint> contactLine() const;
    

};

#endif // PSEUDOPOTEQHANDLER
