#ifndef PDFFIELD_H
#define PDFFIELD_H

#include <latticeField.H>


class pdfField : public latticeField {


protected:


    /** Field values */

    std::vector< std::vector<scalar> > field;


    /** Read field using ensight format */

    const void read();    
    


public:

    
    /** Default constructor */

    pdfField( const latticeMesh& m, timeOptions& t, const std::string& nm, const IO iopt, const IO oopt );


    /** Default destructor */

    ~pdfField();



    /** Access operator */
    
    std::vector<scalar>& operator[] (uint i) {  return field[i];  }

    
    /** Synchronization across procceses */

    const void sync();


    /** Write field using ensight format */

    const void write() const;


};

#endif // PDFFIELD_H