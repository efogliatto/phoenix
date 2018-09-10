#ifndef LBEQUATION_H
#define LBEQUATION_H

#include <scalarField.H>

#include <vectorField.H>

#include <pdfField.H>

#include <map>


class lbEquation {

protected:

    
    /** Equation model name */

    std::string ename;


    /** Reference to external scalar fields */

    std::map<std::string, scalarField*> scfields;


    /** Reference to external vector fields */

    std::map<std::string, vectorField*> vfields;


    /** Reference to external pdf field */

    pdfField* pdffield;



public:


    /** Default constructor */

    lbEquation( const std::string& name );


    /** Default destructor */

    ~lbEquation();



    /** Collision process */

    virtual const void collision() = 0;


    /** Streamming process */

    virtual const void streamming() = 0;


    /** Update macroscopic field */

    virtual const void updateMacroField() = 0;

    

};

#endif // LBEQUATION_H
