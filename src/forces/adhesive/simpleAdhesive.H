#ifndef SIMPLEADHESIVE_H
#define SIMPLEADHESIVE_H

#include <adhesiveForce.H>

class simpleAdhesive : public adhesiveForce {

protected:

    
    /** Interaction weights */

    std::vector<scalar> _weights;

    
    /** Interaction potential */

    const scalar potential( const scalar& rho, const scalar& T, const scalar& cs2 ) const;


    /** Adhesive constant */

    scalar _Gads;


public:

    /** Constructor */

    simpleAdhesive( const std::string& dictName, const std::string& eqName, const latticeMesh& mesh );


    /** Destructor */

    ~simpleAdhesive();


    /** Force at specific node */

    virtual const void force( const uint& i, const scalarField& rho, const scalarField& T,  scalar f[3] ) const;

    
    

};

#endif // SIMPLEADHESIVE_H
