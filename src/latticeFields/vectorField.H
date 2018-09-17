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



    /** Access operator */
    
    std::vector<scalar>& operator[] (uint i) {  return field[i];  }

    const scalar at(const uint& i, const uint& j) const{ return field[i][j]; }

    
    /** Synchronization across procceses */

    const void sync();

    
    /** Start sync */

    const void startSync();


    /** End sync */

    const void endSync();


    
    /** Write field using ensight format */

    const void write() const;        
    

};

#endif // VECTORFIELD_H
