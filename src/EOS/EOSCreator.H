#ifndef EOSCREATOR_H
#define EOSCREATOR_H

#include <vanDerWaals.H>


class EOSCreator {


public:


    // Create with name
    
    EOS* create( const std::string& fname );
    

};


#endif // EOSCREATOR_H
