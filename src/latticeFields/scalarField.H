#ifndef SCALARFIELD_H
#define SCALARFIELD_H

#include <latticeField.H>


class scalarField : public latticeField {


protected:


    /** Field values */

    std::vector<scalar> field;



    /** Read field using ensight format */

    const void read();    
    


public:

    
    /** Default constructor */

    scalarField( const latticeMesh& m, timeOptions& t, const std::string& nm, const IO iopt, const IO oopt );


    /** Default destructor */

    ~scalarField();



    /** Access operator */
    
    scalar& operator[] (uint i) {  return field[i];  }

    
    /** Synchronization across procceses */

    const void sync();


    /** Write field using ensight format */

    const void write() const;    
    
    

};

#endif // SCALARFIELD_H
