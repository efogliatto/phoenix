#ifndef HSCREATOR_H
#define HSCREATOR_H

#include <markusHaziHS.H>

#include <liHS.H>


class hsCreator {

public:

    
    /** Source creator */
    
    heatSource* create( const std::string& dictName, const std::string& eqName, const latticeMesh& mesh, timeOptions& Time );


};

#endif // HSCREATOR_H
