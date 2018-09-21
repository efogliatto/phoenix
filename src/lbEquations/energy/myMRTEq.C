#include <myMRTEq.H>


/** Default constructor */

myMRTEq::myMRTEq( const std::string& name,
		  const latticeMesh& mesh_,
		  timeOptions& Time_,
		  pdfField& pdf_,
		  const scalarField& rho_,
		  const vectorField& U_,
		  scalarField& T_)

    : energyEquation(name, mesh_, Time_, pdf_, rho_, U_, T_) {


    // Read constants

    dictionary dict("properties/macroProperties");

    _a1 = dict.lookUp<scalar>( name + "/HeatSource/Constants/alpha_1");

    _a2 = dict.lookUp<scalar>( name + "/HeatSource/Constants/alpha_2");
    
    

}



/** Default destructor */

myMRTEq::~myMRTEq() {}



/** Collision process */

const void myMRTEq::collision() {}


/** Set pdf to equilibrium values */

const void myMRTEq::setEquilibrium() {}
