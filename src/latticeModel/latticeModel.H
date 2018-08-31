/*
  Abstract base class for DdQq models. Only provides interface
 */


#ifndef LATTICEMODEL_H
#define LATTICEMODEL_H


#include <vector>

#include <string>

#include <dataTypes.H>




class latticeModel {



protected:

  
  /** Dimension */

  uint _d;


  /** Number of velocities */

  uint _q;


  /** Discrete lattice velocities */
  
  std::vector<int> _lvel;


  /** Lattice weights */

  std::vector<scalar> _omega;
  
    
  /** Reverse direction */
  
  std::vector<uint> _reverse;
  


  /** Squared sound speed */
  
  scalar _cs2;
  




  


public:

  
  // Constructors and destructors

  /** Default constructor */

  latticeModel();

  
  /** Default destructor */

  ~latticeModel();




  // Access members
 
  
  /** Dimension */
  
  virtual const uint& d() const = 0;

  
  /** Number of lattice velocities */
  
  virtual const uint q() const = 0;

  
  /** Square of sound speed */
  
  virtual const scalar& cs2() const = 0;


  /** Discrete lattice velocities */
  
  const std::vector<int>& lvel() const;
  
  
  /** Lattice weights */
  
  virtual const std::vector<scalar>& omega() const = 0;
  

  /** Reverse indices */
  
  virtual const std::vector<uint>& reverse() const = 0;
  

};





#endif // LATTICEMODEL_H
