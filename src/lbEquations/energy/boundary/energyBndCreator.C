#include <energyBndCreator.H>

using namespace std;


energyBndCond* energyBndCreator::create(const std::string& eqName,
					const std::string& bdName,
					const latticeMesh& mesh,
					const scalarField& rho,
					const scalarField& T,
					const vectorField& U,
					pdfField& pdf) {

    
    // Load model name from dictionary

    dictionary dict("start/boundaries");

    string btype = dict.lookUp<string>( eqName + "/" + bdName + "/type" );

    
    if( btype == "fixedT" ) {

    	return new energyFixedT( eqName, bdName, mesh, rho, T, U, pdf );

    }

    

    else {


	if( btype == "periodic" ) {

	    return new energyPeriodic( eqName, bdName, mesh, rho, T, U, pdf );

	}

	else {


	    // Default
    
	    cout << endl << " [ERROR]  Boundary condition " << btype << " not available for energy model" << endl << endl;

	    exit(1);

	}
    
    }

    
    return 0;


}
