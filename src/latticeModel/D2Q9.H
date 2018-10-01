#ifndef D2Q9_H
#define D2Q9_H


#include <latticeModel.H>

class D2Q9 : public latticeModel {


    /* ----------------------  Public member functions ----------------------  */

public:

  
    // Constructors and destructors

  
    /** Default constructor */
  
    D2Q9();

  
    /** Default destructor */
  
    ~D2Q9();


    // Acces members

    /** Dimension */
  
    const uint& d() const;


    /** Lattice type */

    const latticeType type() const { return latticeType::D2Q9; }
    
};

#endif // D2Q9_H
