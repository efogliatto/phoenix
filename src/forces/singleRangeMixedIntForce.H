#ifndef SINGLERANGEMIXEDINTFORCE_H
#define SINGLERANGEMIXEDINTFORCE_H

#include <singleRangeIntForce.H>


class singleRangeMixedIntForce : public singleRangeIntForce {


protected:

    
    /** Tunable constant */

    scalar _beta;
    


public:


    /** Constructor */

    singleRangeMixedIntForce( const std::string& dictName, const std::string& eqName, const latticeMesh& mesh, timeOptions& Time );


    /** Destructor */

    ~singleRangeMixedIntForce();


    /** Update force field */

    void update( scalarField& rho, scalarField& T );    
    

};

#endif // SINGLERANGEMIXEDINTFORCE_H
