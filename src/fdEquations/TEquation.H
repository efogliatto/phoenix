#ifndef TEQUATION_H
#define TEQUATION_H

#include <scalarField.H>

#include <vectorField.H>



class TEquation {

protected:

    
    /** Reference to lattice mesh */

    const latticeMesh& _mesh;


    /** Reference to time options */

    timeOptions& _Time;
    
 


public:


    /** Default constructor */

    TEquation( const latticeMesh& mesh, timeOptions& Time, scalarField& T );


    /** Default destructor */

    ~TEquation();



};

#endif // TEQUATION_H
