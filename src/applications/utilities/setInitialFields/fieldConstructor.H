#ifndef FIELDCONSTRUCTOR_H
#define FIELDCONSTRUCTOR_H

#include <dataTypes.H>

#include <scalarField.H>

#include <vectorField.H>

#include <pdfField.H>


template<typename T, typename a>
class fieldConstructor {

public:

    /** Default constructor */

    fieldConstructor();


    /** Default destructor */

    ~fieldConstructor();


    /** Field value at node */
    
    const a field(const T& field, const uint& id) const;

    

};


// Default destructor
template<typename T, typename a>
fieldConstructor<T,a>::fieldConstructor() {}

// Default destructor
template<typename T, typename a>
fieldConstructor<T,a>::~fieldConstructor() {}


#endif // FIELDCONSTRUCTOR_H
