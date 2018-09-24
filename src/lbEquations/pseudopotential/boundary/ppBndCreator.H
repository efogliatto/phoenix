#ifndef PPBNDCREATOR_H
#define PPBNDCREATOR_H

#include <ppFixedU.H>

#include <ppPeriodic.H>

class ppBndCreator {

public:

    ppBndCond* create( const std::string& eqName,
		       const std::string& bdName,
		       const latticeMesh& mesh,
		       const scalarField& rho,
		       const scalarField& T,
		       const vectorField& U,
		       pdfField& pdf );
    
};

#endif // PPBNDCREATOR_H
