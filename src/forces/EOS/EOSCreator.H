#ifndef EOSCREATOR_H
#define EOSCREATOR_H

#include <vanDerWaals.H>

#include <CarnahanStarling.H>

#include <PengRobinson.H>

#include <piecewiseLinear.H>

#include <map>



class EOSCreator {

protected:


    // EOS types

    enum class eosType{ VdW, CS, PR, PWL };

    std::map<std::string, eosType> _eosMapType;

    

public:


    // Create with name
    
    EOS* create( const std::string& dictName, const std::string& eqName );
    

};


#endif // EOSCREATOR_H
