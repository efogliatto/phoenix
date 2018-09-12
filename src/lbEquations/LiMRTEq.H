#ifndef LIMRTEQ_H
#define LIMRTEQ_H

#include <pseudoPotEquation.H>

class LiMRTEq : public pseudoPotEquation {

public:

    
    /** Default constructor */

    LiMRTEq( const std::string& name,
	     const latticeMesh& mesh_,
	     timeOptions& Time_,
	     pdfField pdf_,
	     scalarField& rho_,
	     vectorField& U_,
	     scalarField& T_);
    

    /** Default destructor */

    ~LiMRTEq();



    /** Collision process */

    const void collision();

    
    /** Set pdf to equilibrium values */

    const void setEquilibrium();
    

    

};

#endif // LIMRTEQ_H
