#include <energyPeriodic.H>

/** Constructor */

energyPeriodic::energyPeriodic( const std::string& eqName,
				const std::string& bdName,
				const std::vector<uint>& nodes,
				const scalarField& rho,
				const scalarField& T,
				const vectorField& U,
				pdfField& pdf )
    
    : energyBndCond(nodes, rho, T, U, pdf) {}


/** Destructor */

energyPeriodic::~energyPeriodic() {}


/** Update pdf field */

void energyPeriodic::update( const latticeMesh& mesh ) {}
