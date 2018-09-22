#include <energyPeriodic.H>

/** Constructor */

energyPeriodic::energyPeriodic( const std::string& eqName,
				const std::string& bdName,
				const latticeMesh& mesh,
				const scalarField& rho,
				const scalarField& T,
				const vectorField& U,
				pdfField& pdf )
    
    : energyBndCond(mesh, rho, T, U, pdf, bdName) {}


/** Destructor */

energyPeriodic::~energyPeriodic() {}


/** Update pdf field */

void energyPeriodic::update( const energyEquation* eeq ) {}
