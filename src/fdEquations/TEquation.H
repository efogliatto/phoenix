#ifndef TEQUATION_H
#define TEQUATION_H

#include <scalarField.H>

#include <vectorField.H>

#include <EOSCreator.H>



class TEquation {

protected:

    
    /** Reference to lattice mesh */

    const latticeMesh& _mesh;


    /** Reference to time options */

    timeOptions& _Time;


    /** Specific Heat */

    scalar _Cv;


    /** Diffisivity or conductivity*/

    scalar _lambda;

    
    /** Conductivity model type */

    enum class thmodel { CC, CD };

    thmodel thermal;


    /** Equation of state */

    EOS* eos;
    
 


public:


    /** Default constructor */

    TEquation( const latticeMesh& mesh, timeOptions& Time, scalarField& T );


    /** Default destructor */

    ~TEquation();


    /** Update scalar field using equation info */

    const void rval( scalarField& field, const scalarField& rho, const vectorField& U, const scalarField& T );



};

#endif // TEQUATION_H
