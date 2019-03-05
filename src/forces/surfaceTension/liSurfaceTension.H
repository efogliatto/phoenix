#ifndef LIADHESIVE_H
#define LIADHESIVE_H

#include <simpleAdhesive.H>

class liAdhesive : public simpleAdhesive {

    
public:

    /** Constructor */

    liAdhesive( const std::string& dictName, const std::string& eqName, const latticeMesh& mesh );


    /** Destructor */

    ~liAdhesive();


    /** Force at specific node */

    const void force( const uint& i, const scalarField& rho, const scalarField& T,  scalar f[3] ) const;
    

};

#endif // LIADHESIVE_H
