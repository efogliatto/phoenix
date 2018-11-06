#ifndef D3Q15_H
#define D3Q15_H


#include <latticeModel.H>

class D3Q15 : public latticeModel {


    /* ----------------------  Public member functions ----------------------  */

public:

  
    // Constructors and destructors

  
    /** Default constructor */
  
    D3Q15();

  
    /** Default destructor */
  
    ~D3Q15();


    // Acces members

    /** Dimension */
  
    const uint& d() const;


    /** Lattice type */

    const latticeType type() const { return latticeType::D3Q15; }
    
};

#endif // D3Q15_H
