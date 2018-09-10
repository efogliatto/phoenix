#ifndef VECTORFIELD_H
#define VECTORFIELD_H

#include <latticeField.H>


class vectorField : public latticeField {


protected:


    /** Field values */

    std::vector< std::vector<scalar> > field;



    /** Read field using ensight format */

    const void read();    
    


public:

    
    /** Default constructor */

    vectorField( const latticeMesh& m, timeOptions& t, const std::string& nm, const IO iopt, const IO oopt );


    /** Default destructor */

    ~vectorField();



    // Acces members. Forced interface

    /** Synchronization across procceses */

    const void sync();


    /** Write field using ensight format */

    const void write() const;        
    

};

#endif // VECTORFIELD_H
