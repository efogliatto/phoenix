#ifndef ADSFORCECREATOR_H
#define ADSFORCECREATOR_H

#include <noAds.H>

#include <phiBasedMod.H>


class adsForceCreator {

    
protected:

    /** Force type types */

    enum class adsType{ NONE, PHI_BASED_MOD };

    std::map<std::string, adsType> _adsMapType;

    
    
public:

    /** Force creator */
    
    adhesiveForce* create( const std::string& dictName, const std::string& eqName, const latticeMesh& mesh, timeOptions& Time, const interactionForce* Fi );    

};


#endif // ADSFORCECREATOR_H
