#ifndef LIMRTEQ_H
#define LIMRTEQ_H

#include <pseudoPotEquation.H>

class LiMRTEq : public pseudoPotEquation {

protected:


    /** Sigma constant */

    scalar _sigma;

    
public:

    
    /** Constructor */

    LiMRTEq( const std::string& name,
	     const latticeMesh& mesh_,
	     timeOptions& Time_,
	     pdfField& pdf_,
	     scalarField& rho_,
	     vectorField& U_,
	     scalarField& T_);
    

    /** Default destructor */

    ~LiMRTEq();



    /** Collision process */

    const void collision();

    
    /** Set pdf to equilibrium values */

    const void setEquilibrium();


    /** Compute and update isotropic pressure as scalar field */

    const void pressure( const scalarField& phi, scalarField& p );    
    

    

};

#endif // LIMRTEQ_H
