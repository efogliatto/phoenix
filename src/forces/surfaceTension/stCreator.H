#ifndef STCREATOR_H
#define STCREATOR_H

#include <noSurfaceTension.H>

#include <liSurfaceTension.H>



class stCreator {

public:

    /** Surface Tension creator */
    
    surfaceTension* create( const std::string& dictName, const std::string& eqName, const latticeMesh& mesh );    

};


#endif // STCREATOR_H
