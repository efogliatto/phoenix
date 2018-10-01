/*
  Abstract base class for DdQq models. Only provides interface
 */


#ifndef LATTICEMODEL_H
#define LATTICEMODEL_H


#include <vector>

#include <string>

#include <dataTypes.H>

#include <scalarMatrix.H>




class latticeModel {



protected:

  
    /** Dimension */

    uint _d;


    /** Number of velocities */

    uint _q;


    /** Discrete lattice velocities */
  
    std::vector< std::vector<int> > _lvel;


    /** Lattice weights */

    std::vector<scalar> _omega;
  
    
    /** Reverse direction */
  
    std::vector<uint> _reverse;  


    /** Squared sound speed */
  
    scalar _cs2;


    /** MRT Matrix */

    scalarMatrix M;

    scalarMatrix invM;
    
  




  


public:


    /** Lattice model type */

    enum class latticeType { D2Q9, D3Q15 };

  
    // Constructors and destructors

    /** Default constructor */

    latticeModel();

  
    /** Default destructor */

    ~latticeModel();




    // Access members
 
  
    /** Dimension */
  
    virtual const uint& d() const = 0;


    /** Lattice type */

    virtual const latticeType type() const = 0;

  
    /** Number of lattice velocities */
  
    const uint q() const;

  
    /** Square of sound speed */
  
    const scalar& cs2() const;


    /** Discrete lattice velocities */
  
    const std::vector< std::vector<int> >& lvel() const;
  
  
    /** Lattice weights */
  
    const std::vector<scalar>& omega() const;
  

    /** Reverse indices */
  
    const std::vector<uint>& reverse() const;


    /** Model name */

    const std::string name() const { return "D" + std::to_string(_d) + "Q" + std::to_string(_q);  }


    /** MRT Matrices */

    const scalarMatrix& MRTMatrix() const { return M; }

    const scalarMatrix& MRTInvMatrix() const { return invM; }    

};





#endif // LATTICEMODEL_H
