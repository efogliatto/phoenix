#include <sphere.H>

using namespace std;


// Default constructor

sphere::sphere() : basicShape() {}



// Construct from entry

sphere::sphere( const string& dname, const string& ename ) {

    sphere::readFromEntry(dname, ename );

}



// Default destructor

sphere::~sphere() {}


// Read from entry

void sphere::readFromEntry(const string& dname, const string& ename ) {


    // Load dictionary

    dictionary dict(dname);


    // Read centre

    _centre = dict.lookUp< vector<scalar> >(ename + "/centre");


    // Read radius

    _radius = dict.lookUp<scalar>(ename + "/radius");    
    

}



// Check if point is inside

const bool sphere::isInside(const vector<int>& point) const {

    scalar r0(0);

    for(uint i = 0 ; i < 3 ; i++)
	r0 += (point[i] - _centre[i]) * (point[i] - _centre[i]);

    return r0 <= _radius*_radius;

}



// Bounding box

void sphere::boundingBox( scalar min[3], scalar max[3] ) const {

    for( uint i = 0 ; i < 3 ; i++ ) {
    
	min[i] = _centre[i] - _radius;

	max[i] = _centre[i] + _radius;

    }

}
