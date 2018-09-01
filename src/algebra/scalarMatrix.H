#ifndef SCALARMATRIX_H
#define SCALARMATRIX_H

#include <dataTypes.H>

#include <vector>

#include <scalarVector.H>

#include <iostream>



class scalarMatrix : public std::vector< std::vector<scalar> > {

    
protected:

    /** Cached matrix size */
    
    uint _size;


public:


    /** Default constructor */

    scalarMatrix();



    /** Constructor from vector< vector<scalar> > */

    scalarMatrix( const vector< vector<scalar> >& M );



    /** Default destructor */

    ~scalarMatrix();
    
    

    /** Matrix size */

    const uint& sz() const {   return _size;   }



    /** Update matrix information */

    void update();


    
    /** Matrix - vector multiplication */

    const scalarVector operator* (const scalarVector& V) const;



    /** Overloaded << operator */

    friend std::ostream& operator<<(std::ostream& os, const scalarMatrix& M);
    

};



#endif // SCALARMATRIX_H
