#include <D2Q9.h>

using namespace std;



/* ----------------------  Public member functions ----------------------  */


// Constructors and destructors

// Default constructor
D2Q9::D2Q9() {

    // Model name
    _name = "D2Q9";
    
    // Dimension
    _D = 2;

    // Sound speed
    _cs2 = 1.0 / 3.0;


    // Discrete velocities
    _latticeVel.push_back( Vector3(0,0,0)   );
    _latticeVel.push_back( Vector3(1,0,0)   );
    _latticeVel.push_back( Vector3(0,1,0)   );
    _latticeVel.push_back( Vector3(-1,0,0)  );
    _latticeVel.push_back( Vector3(0,-1,0)  );
    _latticeVel.push_back( Vector3(1,1,0)   );
    _latticeVel.push_back( Vector3(-1,1,0)  );
    _latticeVel.push_back( Vector3(-1,-1,0) );
    _latticeVel.push_back( Vector3(1,-1,0)  );


    // PDF weights
    _omega.push_back( 4.0 / 9.0 );
    _omega.push_back( 1.0 / 9.0);
    _omega.push_back( 1.0 / 9.0 );
    _omega.push_back( 1.0 / 9.0 );
    _omega.push_back( 1.0 / 9.0 );
    _omega.push_back( 1.0 / 36.0 );
    _omega.push_back( 1.0 / 36.0 );
    _omega.push_back( 1.0 / 36.0 );
    _omega.push_back( 1.0 / 36.0 );

    _pdfomega.push_back( 4.0 / 9.0 );
    _pdfomega.push_back( 1.0 / 9.0);
    _pdfomega.push_back( 1.0 / 9.0 );
    _pdfomega.push_back( 1.0 / 9.0 );
    _pdfomega.push_back( 1.0 / 9.0 );
    _pdfomega.push_back( 1.0 / 36.0 );
    _pdfomega.push_back( 1.0 / 36.0 );
    _pdfomega.push_back( 1.0 / 36.0 );
    _pdfomega.push_back( 1.0 / 36.0 );


    // Reverse directions
    _reverse.push_back(0);
    _reverse.push_back(3);
    _reverse.push_back(4);
    _reverse.push_back(1);
    _reverse.push_back(2);
    _reverse.push_back(7);
    _reverse.push_back(8);
    _reverse.push_back(5);
    _reverse.push_back(6);


    // Specular directions
    _specular.push_back(0);
    _specular.push_back(1);
    _specular.push_back(2);
    _specular.push_back(3);
    _specular.push_back(4);
    _specular.push_back(8);
    _specular.push_back(7);
    _specular.push_back(6);
    _specular.push_back(8);


    // MRT matrix
    _M.resize( 9 );
    for(uint i = 0 ; i < 9 ; i++)
	_M[i].resize(9);

    _M[0] = { 1,  1,  1,  1,  1,  1,  1,  1,  1};
    _M[1] = {-4, -1, -1, -1, -1,  2,  2,  2,  2};
    _M[2] = { 4, -2, -2, -2, -2,  1,  1,  1,  1};
    _M[3] = { 0,  1,  0, -1,  0,  1, -1, -1,  1};
    _M[4] = { 0, -2,  0,  2,  0,  1, -1, -1,  1};
    _M[5] = { 0,  0,  1,  0, -1,  1,  1, -1, -1};
    _M[6] = { 0,  0, -2,  0,  2,  1,  1, -1, -1};
    _M[7] = { 0,  1, -1,  1, -1,  0,  0,  0,  0};
    _M[8] = { 0,  0,  0,  0,  0,  1, -1,  1, -1};


    _invM.resize( 9 );
    for(uint i = 0 ; i < 9 ; i++)
	_invM[i].resize(9);

    const double c1 = 1.0/9.0,
	c2 = 1.0 / 6.0,
	c3 = 1.0 / 4.0;

    _invM[0] = {c1,   -c1,      c1,      0,     0,       0,     0,      0,    0};
    _invM[1] = {c1,   -c1/4,   -c1/2,    c2,   -c2,      0,     0,      c3,   0};
    _invM[2] = {c1,   -c1/4,   -c1/2,    0,     0,       c2,   -c2,    -c3,   0};
    _invM[3] = {c1,   -c1/4,   -c1/2,   -c2,    c2,      0,     0,      c3,   0};
    _invM[4] = {c1,   -c1/4,   -c1/2,    0,     0,      -c2,    c2,    -c3,   0};
    _invM[5] = {c1,    c1/2,    c1/4,    c2,    c2/2,    c2,    c2/2,   0,    c3};
    _invM[6] = {c1,    c1/2,    c1/4,   -c2,   -c2/2,    c2,    c2/2,   0,   -c3};
    _invM[7] = {c1,    c1/2,    c1/4,   -c2,   -c2/2,   -c2,   -c2/2,   0,    c3};
    _invM[8] = {c1,    c1/2,    c1/4,    c2,    c2/2,   -c2,   -c2/2,   0,   -c3};




    _principal.resize( 3 );
    for(uint i = 0 ; i < 3 ; i++)
    	_principal[i].resize(2);

    _principal[0][0] = 3;
    _principal[0][1] = 1;
    _principal[1][0] = 4;
    _principal[1][1] = 2;
    _principal[2][0] = 0;
    _principal[2][1] = 0;



    // Pseudopotential weights

    _weights.push_back( 0 );
    _weights.push_back(1.0/3);
    _weights.push_back(1.0/3);
    _weights.push_back(1.0/3);
    _weights.push_back(1.0/3);
    _weights.push_back(1.0/12);
    _weights.push_back(1.0/12);
    _weights.push_back(1.0/12);
    _weights.push_back(1.0/12);
    
}


// Default destructor
D2Q9::~D2Q9() {}



// Acces members
const uint& D2Q9::D() const {
    return _D;
}



// Main index
const bool D2Q9::is_principal(const uint& id) const {
    return ( (id <= 4)  &&  (id > 0) ) ? true : false;    
}
