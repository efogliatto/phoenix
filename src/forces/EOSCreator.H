#ifndef EOSCREATOR_H
#define EOSCREATOR_H

#include <vanDerWaals.H>


class EOSCreator {


public:


    // Create with name
    
    EOS* create( const std::string& dictName, const std::string& eqName );
    

};


#endif // EOSCREATOR_H
