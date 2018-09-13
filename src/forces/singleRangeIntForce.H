#ifndef SINGLERANGEINTFORCE_H
#define SINGLERANGEINTFORCE_H

#include <interactionForce.H>


class singleRangeIntForce : public interactionForce {


protected:


    /** Interaction weights */

    std::vector<scalar> _weights;


public:


    /** Constructor */

    singleRangeIntForce( const std::string& dictName, const std::string& eqName, const latticeMesh& mesh, timeOptions& Time );


    /** Destructor */

    ~singleRangeIntForce();


    /** Update force field */

    void update( scalarField& rho, scalarField& T );    
    

};

#endif // SINGLERANGEINTFORCE_H
