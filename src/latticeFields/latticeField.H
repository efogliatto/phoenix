#ifndef LATTICEFIELD_H
#define LATTICEFIELD_H

#include <latticeMesh.H>

#include <timeOptions.H>


class latticeField {


protected:


    /** Field name */

    const std::string name;
    
    
    /** Reference to lattice mesh */

    const latticeMesh& mesh;


    /** Reference to time options */

    timeOptions& Time;




    /** Read field using ensight format */

    virtual const void read() = 0;    
    


public:

    /** Default constructor */

    latticeField( const latticeMesh& m, timeOptions& t, const std::string& nm );


    /** Default destructor */

    ~latticeField();



    // Acces members. Forced interface

    /** Synchronization across procceses */

    virtual const void sync() = 0;


    /** Write field using ensight format */

    virtual const void write() const= 0;    
    
    

};

#endif // LATTICEFIELD_H
