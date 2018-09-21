#ifndef MYMRTEQ_H
#define MYMRTEQ_H

#include <energyEquation.H>

#include <hsCreator.H>


class myMRTEq : public energyEquation {


protected:

    
    /** Constants */

    scalar _a1;

    scalar _a2;
  



public:

    
    /** Default constructor */

    myMRTEq( const std::string& name,
	     const latticeMesh& mesh_,
	     timeOptions& Time_,
	     pdfField& pdf_,
	     const scalarField& rho_,
	     const vectorField& U_,
	     scalarField& T_);


    /** Default destructor */

    ~myMRTEq();



    /** Collision process */

    const void collision();


    /** Set pdf to equilibrium values */

    const void setEquilibrium();

    
};


#endif // MYMRTEQ_H
