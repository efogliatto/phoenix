#ifndef PPBNDCREATOR_H
#define PPBNDCREATOR_H

#include <ppFixedU.H>

#include <ppPeriodic.H>

class ppBndCreator {

public:

    ppBndCond* create( const std::string& eqName, const std::string& bdName, const std::vector<uint>& nodes );

};

#endif // PPBNDCREATOR_H
