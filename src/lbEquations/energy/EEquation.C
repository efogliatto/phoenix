#include <EEquation.H>

using namespace std;


energyEquation* EEquation::create( const std::string& name,
				   const latticeMesh& mesh_,
				   timeOptions& Time_,
				   pdfField& pdf_,
				   const scalarField& rho_,
				   const vectorField& U_,
				   scalarField& T_) {

    
    // Load model name from dictionary

    dictionary dict("properties/macroProperties");

    string etype = dict.lookUp<string>( name + "/LBModel/type" );

    
    if( etype == "myMRT" ) {

    	return new myMRTEq(name, mesh_, Time_, pdf_, rho_, U_, T_);

    }

    else {


	if( etype == "LiEnergyMRT" ) {

	    return new LiEnergyMRTEq(name, mesh_, Time_, pdf_, rho_, U_, T_);

	}

	else {

	
	    // Default
    
	    cout << endl << " [ERROR]  LB equation type " << etype << " not available as energy model" << endl << endl;

	    exit(1);

	}
    
    }

    
    return 0;


}
