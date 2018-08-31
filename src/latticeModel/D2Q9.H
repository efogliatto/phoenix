#ifndef D2Q9_H
#define D2Q9_H


#include <basicLBModel.h>

class D2Q9 : public basicLBModel {


    /* ----------------------  Public member functions ----------------------  */

public:

    // Constructors and destructors

    // Default constructor
    D2Q9();

    // Default destructor
    ~D2Q9();


    // Acces members
    const uint& D() const;

    // Main index
    const bool is_principal(const uint& id) const;
    
};

#endif // D2Q9_H
