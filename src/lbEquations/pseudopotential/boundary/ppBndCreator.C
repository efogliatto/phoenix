#include <ppBndCreator.H>

using namespace std;


ppBndCond* ppBndCreator::create( const string& eqName,
				 const string& bdName,
				 const latticeMesh& mesh,
				 const scalarField& rho,
				 const scalarField& T,
				 const vectorField& U,
				 pdfField& pdf) {

    
    // Load model name from dictionary

    dictionary dict("start/boundaries");

    string btype = dict.lookUp<string>( eqName + "/" + bdName + "/type" );

    
    if( btype == "fixedU" ) {

    	return new ppFixedU( eqName, bdName, mesh, rho, T, U, pdf );

    }

    

    else {


	if( btype == "periodic" ) {

	    return new ppPeriodic( eqName, bdName, mesh, rho, T, U, pdf );

	}

	else {

	    if( btype == "outflow" ) {

		return new ppOutflow( eqName, bdName, mesh, rho, T, U, pdf );

	    }

	    else {

		// Default
    
		cout << endl << " [ERROR]  Boundary condition " << btype << " not available for pseudopotential model" << endl << endl;

		exit(1);

	    }

	}
    
    }

    
    return 0;


}