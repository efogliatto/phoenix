#ifndef LATTICEMODELCREATOR_H
#define LATTICEMODELCREATOR_H

#include <D2Q9.H>


class latticeModelCreator {


public:

    latticeModel* create( const std::string& modelName );

};


#endif // LATTICEMODELCREATOR_H
