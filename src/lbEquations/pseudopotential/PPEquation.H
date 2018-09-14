#ifndef PPEQUATION_H
#define PPEQUATION_H

#include <LiMRTEq.H>

class PPEquation {

public:

    pseudoPotEquation* create( const std::string& name,
			       const latticeMesh& mesh_,
			       timeOptions& Time_,
			       pdfField& pdf_,
			       scalarField& rho_,
			       vectorField& U_,
			       scalarField& T_ );

};

#endif // PPEQUATION_H
