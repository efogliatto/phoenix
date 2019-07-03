#ifndef SIMPLIFIEDTEQ_H
#define SIMPLIFIEDTEQ_H

#include <scalarField.H>

#include <vectorField.H>

#include <EOSCreator.H>

#include <algebra.H>



class simplifiedTEq {

protected:

    
    /** Reference to lattice mesh */

    const latticeMesh& _mesh;


    /** Reference to time options */

    timeOptions& _Time;


    /** Specific Heat */

    scalar _Cv;


    /** Model constants */

    scalar alpha_1;

    scalar alpha_2;

    
    /** Energy equation relaxation coefficients */

    std::vector<scalar> _Tau;


    /** Auxiliary coefficient matrix (Qaux = (0.5-Q^{-1})) */

    sparseScalarMatrix Qaux;
    

    /** Equation of state */

    EOS* eos;
    
 


public:


    /** Default constructor */

    simplifiedTEq( const latticeMesh& mesh, timeOptions& Time, scalarField& T );


    /** Default destructor */

    ~simplifiedTEq();


    // /** Update scalar field using equation info */

    // const void rval( scalarField& field, const scalarField& rho, const vectorField& U, const scalarField& T );


    /** Source value at specific node */

    const scalar heatSource( const uint id ) const;



};

#endif // SIMPLIFIEDTEQ_H
