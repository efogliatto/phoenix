#ifndef LBMODELCREATOR_H
#define LBMODELCREATOR_H

#include <D2Q9.h>
#include <D2Q5.h>
#include <D2Q4.h>
#include <D3Q7.h>

class LBModelCreator {


    /* ----------------------  Public member functions ----------------------  */

public:

    basicLBModel* create(const std::string& modelName);

};


#endif // LBMODELCREATOR_H
