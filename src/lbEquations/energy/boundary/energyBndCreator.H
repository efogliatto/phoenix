#ifndef ENERGYBNDCREATOR_H
#define ENERGYBNDCREATOR_H

#include <energyFixedT.H>

#include <energyPeriodic.H>


class energyBndCreator {

public:

    energyBndCond* create(const std::string& eqName,
			  const std::string& bdName,
			  const latticeMesh& mesh,
			  const scalarField& rho,
			  const scalarField& T,
			  const vectorField& U,
			  pdfField& pdf);

};

#endif // ENERGYBNDCREATOR_H
