#ifndef PPBNDCREATOR_H
#define PPBNDCREATOR_H

#include <ppFixedU.H>

#include <ppNEBB.H>

#include <ppGeneralNEBB.H>

#include <ppNEExt.H>

#include <ppPeriodic.H>

#include <ppOutflow.H>

#include <ppOutflowWithNEBB.H>

class ppBndCreator {

    
protected:

    
    // Boundary types

    enum class bdType{ fixedU, periodic, outflow, NEBB, NEExt, G_NEBB, OutflowNEBB };

    std::map<std::string, bdType> _bdMapType;


    
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
