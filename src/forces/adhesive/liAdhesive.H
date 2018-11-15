#ifndef LIADHESIVE_H
#define LIADHESIVE_H

#include <adhesiveForce.H>

class liAdhesive : public adhesiveForce {

protected:

    /** Adhesion coefficient for each boundary */

    std::map<std::string, scalar> _Gads;

    
    /** Interaction weights */

    std::vector<scalar> _weights;


    
    
    /** Interaction potential */

    const scalar potential( const scalar& rho, const scalar& T, const scalar& cs2 ) const;
    


public:

    /** Constructor */

    liAdhesive( const std::string& dictName, const std::string& eqName, const latticeMesh& mesh );


    /** Destructor */

    ~liAdhesive();


    /** Force at specific node */

    const void force( const uint& i, const scalar& rho, const scalar& T,  scalar f[3] ) const;
    

};

#endif // LIADHESIVE_H
