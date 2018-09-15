#include <ppPeriodic.H>

using namespace std;


/** Constructor */

ppPeriodic::ppPeriodic( const std::vector<uint>& nodes ) : ppBndCond(nodes) {}


/** Destructor */

ppPeriodic::ppPeriodic() {}


/** Update pdf field */

void ppPeriodic::update( const latticeMesh& mesh, pdfField& pdf, const scalarField& rho, const scalarField& T, const vectorField& U, const pseudoPotForce* F ) {}
