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

    const scalar& at(uint i) const {  return field[i];  }


    /** Gradient at node */

    const void grad(scalar g[3], const uint& id, const bool inverse = false) const;

    
    /** Synchronization across procceses */

    const void sync();


    /** Start sync */

    const void startSync();


    /** End sync */

    const void endSync();
    

    /** Write field using ensight format */

    const void write() const;


    /** Field average */

    const scalar average() const;


    /** Update field values from time index (not time value) */
    
    const void update( const uint& t );

};

#endif // SCALARFIELD_H
