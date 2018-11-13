#ifndef ENERGYBNDCREATOR_H
#define ENERGYBNDCREATOR_H

#include <energyFixedT.H>

#include <energyNormalHeatFlux.H>

#include <energyPeriodic.H>

#include <energyOutflow.H>

#include <energyFixedTSpot.H>


class energyBndCreator {

protected:

    
    // Boundary types

    enum class bdType{ fixedT, periodic, fixedGradT, normalHeatFlux, outflow, fixedTSpot };

    std::map<std::string, bdType> _bdMapType;
    
    
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
