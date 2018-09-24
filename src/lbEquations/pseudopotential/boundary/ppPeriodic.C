#include <ppPeriodic.H>

using namespace std;


/** Constructor */

ppPeriodic::ppPeriodic( const std::string& eqName,
			const std::string& bdName,
			const latticeMesh& mesh,
			const scalarField& rho,
			const scalarField& T,
			const vectorField& U,
			pdfField& pdf)

    : ppBndCond(mesh, rho, T, U, pdf, bdName) {}


/** Destructor */

ppPeriodic::~ppPeriodic() {}


/** Update pdf field */

void ppPeriodic::update( const pseudoPotEquation* ppeq ) {}
