#ifndef EXTERNALFORCE_H
#define EXTERNALFORCE_H

#include <scalarField.H>

#include <vectorField.H>


class externalForce {

protected:
    

    /** Force magnitude */
    
    std::vector<scalar> _force;


public:

    /** Constructor */

    externalForce( const std::string& dictName, const std::string& eqName );


    /** Destructor */

    ~externalForce();


    /** Force */

    const std::vector<scalar>& force() const { return _force; }


    /** Overload access operator */

    const scalar& operator[] (const uint& i) const { return _force[i]; }

};

#endif // EXTERNALFORCE_H
