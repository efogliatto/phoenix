#ifndef SPARSESCALARMATRIX_H
#define SPARSESCALARMATRIX_H

#include <dataTypes.H>

#include <vector>

#include <iostream>



class sparseScalarMatrix {

    
protected:


    /** Matrix coefficients */

    std::vector<scalar> _values;


    /** Matrix indices */

    std::vector< std::vector<uint> > _idx;

    
    /** Cached matrix size */
    
    uint _size;


public:


    /** Default constructor */

    sparseScalarMatrix();



    /** Constructor from diagonal */

    sparseScalarMatrix( const std::vector<scalar>& V );



    /** Default destructor */

    ~sparseScalarMatrix();
    


    /** Matrix - vector multiplication */

    const void matDotVec (const std::vector<scalar>& V, std::vector<scalar>& res) const;


    
    /** Add element */

    const void addElement( const scalar& val, const uint& i, const uint& j );
    
    

};



#endif // SPARSESCALARMATRIX_H
