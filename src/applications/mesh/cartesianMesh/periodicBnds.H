#ifndef PERIODICBOUNDARIES_H
#define PERIODICBOUNDARIES_H

#include <string>

#include <vector>

#include <dataTypes.H>

struct periodicBnds {

    std::string bd1;

    std::string bd2;

    std::vector<scalar> direction;

};

#endif // PERIODICBOUNDARIES_H
