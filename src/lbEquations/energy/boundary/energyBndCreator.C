#include <energyBndCreator.H>

using namespace std;


energyBndCond* energyBndCreator::create(const std::string& eqName,
					const std::string& bdName,
					const latticeMesh& mesh,
					const scalarField& rho,
					const scalarField& T,
					const vectorField& U,
					pdfField& pdf) {



    // Initialize types

    _bdMapType["fixedT"]              = bdType::fixedT;
    _bdMapType["fixedGradT"]          = bdType::fixedGradT;
    _bdMapType["periodic"]            = bdType::periodic;
    _bdMapType["normalHeatFlux"]      = bdType::normalHeatFlux;
    _bdMapType["normalHeatFluxSpots"] = bdType::normalHeatFluxSpots;
    _bdMapType["outflow"]             = bdType::outflow;
    _bdMapType["fixedTSpots"]         = bdType::fixedTSpots;
    _bdMapType["fixedCosT"]           = bdType::fixedCosT;
    _bdMapType["InamuroFixedT"]       = bdType::InamuroFT;
    _bdMapType["eqFixedT"     ]       = bdType::eqFT;
    _bdMapType["NEExtFixedT"]         = bdType::NEExtFT;
    _bdMapType["eqFixedTSpots"]       = bdType::eqFTSpot;        
    _bdMapType["NEExtFixedTSpots"]    = bdType::NEExtFTSpot;
    
    
    
    // Load model name from dictionary

    dictionary dict("start/boundaries");

    string btype = dict.lookUp<string>( eqName + "/" + bdName + "/type" );



    // Assign boundary condition
    
    if( _bdMapType.find(btype) != _bdMapType.end() ) {
	

	switch( _bdMapType[btype] ) {
	

	case bdType::fixedT:

	    return new energyFixedT( eqName, bdName, mesh, rho, T, U, pdf );

	    break;


	case bdType::periodic:

	    return new energyPeriodic( eqName, bdName, mesh, rho, T, U, pdf );

	    break;


	case bdType::fixedGradT:

	    return new energyFixedGradT( eqName, bdName, mesh, rho, T, U, pdf ); 

	    break;


	case bdType::normalHeatFlux:

	    return new energyNormalHeatFlux( eqName, bdName, mesh, rho, T, U, pdf );

	    break;


	case bdType::normalHeatFluxSpots:

	    return new energyNormalHeatFluxSpot( eqName, bdName, mesh, rho, T, U, pdf );

	    break;
	    
	    
	case bdType::outflow:

	    return new energyOutflow( eqName, bdName, mesh, rho, T, U, pdf );

	    break;


	case bdType::fixedTSpots:

	    return new energyFixedTSpot( eqName, bdName, mesh, rho, T, U, pdf );

	    break;


	case bdType::fixedCosT:

	    return new energyFixedCosT( eqName, bdName, mesh, rho, T, U, pdf );

	    break;


	case bdType::InamuroFT:

	    return new InamuroFixedT( eqName, bdName, mesh, rho, T, U, pdf );

	    break;		    


	case bdType::eqFT:

	    return new eqFixedT( eqName, bdName, mesh, rho, T, U, pdf );

	    break;		    


	case bdType::NEExtFT:

	    return new NEExtFixedT( eqName, bdName, mesh, rho, T, U, pdf );

	    break;



	case bdType::eqFTSpot:

	    return new eqFixedTSpot( eqName, bdName, mesh, rho, T, U, pdf );

	    break;	    
	    

	case bdType::NEExtFTSpot:

	    return new NEExtFixedTSpot( eqName, bdName, mesh, rho, T, U, pdf );

	    break;	    
	    
	    
	}

    
    }

    
    else {
    
	cout << endl << " [ERROR]  Boundary condition " << btype << " not available for energy model" << endl << endl;

	exit(1);

    }    
    
    

    
    return 0;


}
