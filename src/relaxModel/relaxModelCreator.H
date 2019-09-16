#ifndef RELAXMODELCREATOR_H
#define RELAXMODELCREATOR_H

#include <uniformTau.H>

#include <rhoPieceWise.H>

#include <rhoPieceWiseLinear.H>

#include <map>

#include <memory>


class relaxModelCreator {


protected:

    /** Model types */

    enum class relaxType{utau, rhopw, rhopwl};

    std::map<std::string, relaxType> _relaxMapType;


public:

    /** Model creation */
    
    std::unique_ptr<relaxModel> create( const std::string& entry );

};

#endif // RELAXMODELCREATOR_H
