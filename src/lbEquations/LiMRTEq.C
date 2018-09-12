#include <LiMRTEq.H>

using namespace std;


/** Default constructor */

LiMRTEq::LiMRTEq( const string& name,
		  const latticeMesh& mesh_,
		  timeOptions& Time_,
		  pdfField pdf_,
		  scalarField& rho_,
		  vectorField& U_,
		  scalarField& T_) : pseudoPotEquation(name,
						       mesh_,
						       Time_,
						       pdf_,
						       rho_,
						       U_,
						       T_) {}
    

/** Default destructor */

LiMRTEq::~LiMRTEq() {}



/** Collision process */

const void LiMRTEq::collision() {}

    
/** Set pdf to equilibrium values */

const void LiMRTEq::setEquilibrium() {}
