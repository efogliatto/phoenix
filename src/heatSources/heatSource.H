#ifndef HEATSOURCE_H
#define HEATSOURCE_H

#include <vectorField.H>

#include <scalarField.H>

#include <EOSCreator.H>



class heatSource {


protected:


    /** Reference to lattice mesh */

    const latticeMesh& _mesh;

    
    /** Reference to time */

    const timeOptions& _Time;
    
    
    /** Force field */

    scalarField _source;




    


public:


    /** Constructor */

    heatSource( const std::string& dictName, const std::string& eqName, const latticeMesh& mesh, timeOptions& Time );


    /** Destructor */

    ~heatSource();


    /** Source at specific node */

    const scalar source( const uint& i ) { return _source[i]; }
    

    /** Update source field */

    virtual void update( const scalarField& rho, const scalarField& T, const vectorField& U ) = 0;


    /** Synv source field */

    const void sync() { _source.sync(); }
    

};

#endif // HEATSOURCE_H
