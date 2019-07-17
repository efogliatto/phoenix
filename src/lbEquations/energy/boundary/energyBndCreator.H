#ifndef ENERGYBNDCREATOR_H
#define ENERGYBNDCREATOR_H

#include <energyFixedT.H>

#include <InamuroFixedT.H>

#include <NEExtFixedT.H>

#include <eqFixedT.H>

#include <energyNormalHeatFlux.H>

#include <energyPeriodic.H>

#include <energyOutflow.H>

#include <energyFixedTSpot.H>

#include <energyFixedCosT.H>


class energyBndCreator {

protected:

    
    // Boundary types

    enum class bdType{ fixedT, periodic, fixedGradT, normalHeatFlux, outflow, fixedTSpots, fixedCosT, InamuroFT, eqFT, NEExtFT };

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
