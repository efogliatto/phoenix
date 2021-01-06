#ifndef INITIALSCALARFIELD_H
#define INITIALSCALARFIELD_H

#include <scalarField.H>

#include <map>


class initialScalarField {

    
protected:

    
    // Initial types

    enum class ishape{ uniform, box, random, bgsphere, linear, tg };

    std::map<std::string, ishape> _spMapType;
    
    
public:

    
    // Default constructor

    initialScalarField();

    
    // Update field
    
    void updateField( scalarField& field, const latticeMesh& mesh, const std::string& fname, const std::string& sptype );

};

#endif // INITIALSCALARFIELD_H
