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

    std::vector<uint> _idx0;

    std::vector<uint> _idx1;

    
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


    /** Access to elements */

    const std::vector<scalar>& values() const {return _values;}
    const std::vector<uint>& idx0() const {return _idx0;}
    const std::vector<uint>& idx1() const {return _idx1;}        
    
    

};



#endif // SPARSESCALARMATRIX_H
