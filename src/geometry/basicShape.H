#ifndef BASICSHAPE_H
#define BASICSHAPE_H

#include <dataTypes.H>

#include <dictionary.H>


class basicShape {

public:

    // Default constructor

    basicShape();


    // Construct from entry

    basicShape( const std::string& dname, const std::string& ename );
    

    // Default destructor

    ~basicShape();


    // Read from entry

    virtual void readFromEntry(const std::string& dname, const std::string& ename ) = 0;


    // Check if point is inside

    virtual const bool isInside(const std::vector<int>& point) const = 0;


    // Bounding box

    virtual void boundingBox( scalar min[3], scalar max[3] ) const = 0;

};

#endif // BASICSHAPE_H
