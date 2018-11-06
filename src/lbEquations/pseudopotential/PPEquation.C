#include <PPEquation.H>

using namespace std;

pseudoPotEquation* PPEquation::create( const std::string& name,
				       const latticeMesh& mesh_,
				       timeOptions& Time_,
				       pdfField& pdf_,
				       scalarField& rho_,
				       vectorField& U_,
				       scalarField& T_ ) {

    
    // Load model name from dictionary

    dictionary dict("properties/macroProperties");

    string etype = dict.lookUp<string>( name + "/LBModel/type" );

    
    if( etype == "LiMRT" ) {

    	return new LiMRTEq(name, mesh_, Time_, pdf_, rho_, U_, T_);

    }

    else {

	if( etype == "XuMRT" ) {

	    return new XuMRTEq(name, mesh_, Time_, pdf_, rho_, U_, T_);

	}

	else {
	

	    // Default
    
	    cout << endl << " [ERROR]  LB equation type " << etype << " not available as pseudopotential model" << endl << endl;

	    exit(1);

	}
    
    }

    
    return 0;

}
