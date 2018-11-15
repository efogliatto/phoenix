#ifndef ADSFORCECREATOR_H
#define ADSFORCECREATOR_H

#include <noAds.H>

#include <liAdhesive.H>


class adsForceCreator {

public:

    /** Force creator */
    
    adhesiveForce* create( const std::string& dictName, const std::string& eqName, const latticeMesh& mesh );    

};


#endif // ADSFORCECREATOR_H
