#ifndef PPBNDCREATOR_H
#define PPBNDCREATOR_H

#include <ppFixedU.H>

#include <ppNEBB.H>

#include <ppPeriodic.H>

#include <ppOutflow.H>

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
