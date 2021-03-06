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


    /** The copy-constructor */

    scalarMatrix( const scalarMatrix& M );
    
    

    /** Matrix size */

    const uint& sz() const {   return _size;   }



    /** Update matrix information */

    void update();
   


    /** Overloaded << operator */

    friend std::ostream& operator<<(std::ostream& os, const scalarMatrix& M);


    /** Overloaded copy operator */

    scalarMatrix& operator= (const scalarMatrix& M);


    /** Set from array */

    void setFromArray( const std::vector< std::vector<scalar> >& M );



    /** Matrix - vector multiplication */

    const void matDotVec (const scalarVector& V, scalarVector& res) const;


    /** Matrix - vector multiplication */

    const void matDotVec (const std::vector<scalar>& V, std::vector<scalar>& res) const;
    
    

};



#endif // SCALARMATRIX_H
