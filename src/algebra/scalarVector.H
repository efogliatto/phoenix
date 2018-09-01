#ifndef SCALARVECTOR_H
#define SCALARVECTOR_H

#include <dataTypes.H>

#include <vector>

#include <iostream>

// #include <scalarMatrix.H>


class scalarVector : public std::vector<scalar> {

    
protected:

    /** Cached vector size */
    
    uint _size;


public:



    /** Default constructor */

    scalarVector();



    /** Constructor from vector< vector<scalar> > */

    scalarVector( const vector<scalar>& V );



    /** Constructor with size and default value */

    scalarVector(uint sz, uint df);    
    


    /** Default destructor */

    ~scalarVector();

    

    /** Vector size */

    const uint& sz() const {   return _size;   }



    /** Update vector information */

    void update();


    
    /** Overloaded << operator */

    friend std::ostream& operator<<(std::ostream& os, const scalarVector& V);
    

    
    // /** Matrix - vector multiplication */
    
    // friend scalarVector& operator*(const scalarMatrix& M, const scalarVector& V);
    
};





#endif // SCALARVECTOR_H
