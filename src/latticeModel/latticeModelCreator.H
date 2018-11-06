#ifndef LATTICEMODELCREATOR_H
#define LATTICEMODELCREATOR_H

#include <D2Q9.H>

#include <D3Q15.H>


class latticeModelCreator {


public:


    // Create with name
    
    latticeModel* create( const std::string& modelName );


    // Create with d and q

    latticeModel* create( const uint& d, const uint& q );
    

};


#endif // LATTICEMODELCREATOR_H
