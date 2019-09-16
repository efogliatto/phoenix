#ifndef RELAXMODEL_H
#define RELAXMODEL_H

#include <string>

#include <dataTypes.H>


class relaxModel {

protected:

    /** Model name */

    std::string _name;


public:

    /** Default constructor */

    relaxModel( const std::string& entry );


    /** Destructor */

    ~relaxModel();


    /** Density-dependent relaxation factor */

    virtual const scalar tau( const scalar rho, const uint i = 0 );
    

};


#endif // RELAXMODEL_H
