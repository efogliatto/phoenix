#ifndef SPHERE_H
#define SPHERE_H

#include <basicShape.H>

class sphere : public basicShape {

    
protected:

    // Sphere centre

    std::vector<scalar> _centre;


    // Radius

    scalar _radius;


public:

    
    // Default constructor

    sphere();


    // Construct from entry

    sphere( const std::string& dname, const std::string& ename );
    

    // Default destructor

    ~sphere();


    // Read from entry

    void readFromEntry( const std::string& dname, const std::string& ename );


    // Check if point is inside

    const bool isInside(const std::vector<int>& point) const;


    // Bounding box

    void boundingBox( scalar min[3], scalar max[3] ) const;


    // Centre

    const std::vector<scalar>& centre() const {return _centre;}


    // Radius

    const scalar radius() const {return _radius;} 

};

#endif // SPHERE_H
