#ifndef EEQUATION_H
#define EEQUATION_H

#include <myMRTEq.H>

#include <GongSRTEq.H>

#include <LiEnergyMRTEq.H>

class EEquation {

public:

    energyEquation* create( const std::string& name,
			    const latticeMesh& mesh_,
			    timeOptions& Time_,
			    pdfField& pdf_,
			    const scalarField& rho_,
			    const vectorField& U_,
			    scalarField& T_);

};

#endif // EEQUATION_H
