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


    /** New temperature */

    scalarField Tnew;
    
 


public:


    /** Default constructor */

    simplifiedTEq( const latticeMesh& mesh, timeOptions& Time, scalarField& T );


    /** Default destructor */

    ~simplifiedTEq();

    
    /** Predictor step (convection + source) */

    const void predictor( scalarField& Tstar, const scalarField& rho, const vectorField& U, const scalarField& T );


    /** Corrector step (diffusion) */

    const void corrector( scalarField& T, const scalarField& rho, const vectorField& U, const scalarField& Tstar );
    

    /** Source value at specific node */

    const scalar heatSource( const scalarField& rho, const scalarField& T, const vectorField& U, const uint id ) const;


    /** Equilibrium value at specific node and velocity */

    const scalar neq( const scalarField& T, const vectorField& U, const uint id, const uint vid );

    



};

#endif // SIMPLIFIEDTEQ_H
