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


    /** MPI send buffer */

    scalar** sbuf;


    /** MPI recv buffer */

    scalar** rbuf;


    /** MPI requests */

    MPI_Request _request[100];

    MPI_Status _status[100];    
    
    int _nreq;

    


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


    /** Start sync */

    virtual const void startSync() = 0;


    /** End sync */

    virtual const void endSync() = 0;


    /** Write field using ensight format */

    virtual const void write() const= 0;


    /** Field name */

    const std::string fieldName() const { return name; }
    

    /** Reference to time */

    const timeOptions& time() const {return Time;}
    

};

#endif // LATTICEFIELD_H
