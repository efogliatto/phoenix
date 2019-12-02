#include <ppBndCreator.H>

using namespace std;


ppBndCond* ppBndCreator::create( const string& eqName,
				 const string& bdName,
				 const latticeMesh& mesh,
				 const scalarField& rho,
				 const scalarField& T,
				 const vectorField& U,
				 pdfField& pdf) {



    // Initialize types

    _bdMapType["fixedU"]   = bdType::fixedU;
    _bdMapType["periodic"] = bdType::periodic;
    _bdMapType["outflow"]  = bdType::outflow;
    _bdMapType["NEBB"]     = bdType::NEBB;
    _bdMapType["NEExt"]    = bdType::NEExt;
    _bdMapType["generalNEBB"]   = bdType::G_NEBB;
    _bdMapType["outflowWithNEBB"]   = bdType::OutflowNEBB;	

    
    
    // Load model name from dictionary

    dictionary dict("start/boundaries");

    string btype = dict.lookUp<string>( eqName + "/" + bdName + "/type" );



    // Assign boundary condition
    
    if( _bdMapType.find(btype) != _bdMapType.end() ) {
	

	switch( _bdMapType[btype] ) {

	case bdType::fixedU:

	    return new ppFixedU( eqName, bdName, mesh, rho, T, U, pdf );

	    break;

	    
	case bdType::periodic:

	    return new ppPeriodic( eqName, bdName, mesh, rho, T, U, pdf );

	    break;


	case bdType::outflow:

	    new ppOutflow( eqName, bdName, mesh, rho, T, U, pdf );

	    break;


	case bdType::NEBB:

	    return new ppNEBB( eqName, bdName, mesh, rho, T, U, pdf );

	    break;


	case bdType::NEExt:

	    new ppNEExt( eqName, bdName, mesh, rho, T, U, pdf );

	    break;


	case bdType::G_NEBB:

	    return new ppGeneralNEBB( eqName, bdName, mesh, rho, T, U, pdf );

	    break;

	    
	case bdType::OutflowNEBB:

	    return new ppOutflowWithNEBB( eqName, bdName, mesh, rho, T, U, pdf );

	    break;	    
	    
	    
	}

    }


    else {
    
	cout << endl << " [ERROR]  Boundary condition " << btype << " not available for pseudopotential model" << endl << endl;

	exit(1);

    }
    

    
    // if( btype == "fixedU" ) {

    // 	return new ppFixedU( eqName, bdName, mesh, rho, T, U, pdf );

    // }

    

    // else {


    // 	if( btype == "periodic" ) {

    // 	    return new ppPeriodic( eqName, bdName, mesh, rho, T, U, pdf );

    // 	}

    // 	else {

    // 	    if( btype == "outflow" ) {

    // 		return new ppOutflow( eqName, bdName, mesh, rho, T, U, pdf );

    // 	    }

    // 	    else {

    // 		if( btype == "NEBB" ) {

    // 		    return new ppNEBB( eqName, bdName, mesh, rho, T, U, pdf );

    // 		}

    // 		else {
	
    // 		    if( btype == "NEExt" ) {

    // 			return new ppNEExt( eqName, bdName, mesh, rho, T, U, pdf );

    // 		    }

    // 		    else {

    // 			// Default
    
    // 			cout << endl << " [ERROR]  Boundary condition " << btype << " not available for pseudopotential model" << endl << endl;

    // 			exit(1);

    // 		    }

    // 		}

    // 	    }

    // 	}
    
    // }

    
    return 0;


}
