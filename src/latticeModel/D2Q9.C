#include <D2Q9.H>

using namespace std;



/* ----------------------  Public member functions ----------------------  */


// Constructors and destructors

/** Default constructor */

D2Q9::D2Q9() {

    
  // Dimension

  _d = 2;


  // Velocities

  _q = 9;
  

  // Sound speed
  
  _cs2 = 1.0 / 3.0;


  // Discrete velocities
  
  _lvel ={  {  0,  0,  0 },
	    {  1,  0,  0 },
	    {  0,  1,  0 },
	    { -1,  0,  0 },
	    {  0, -1,  0 },
	    {  1,  1,  0 },
	    { -1,  1,  0 },
	    { -1, -1,  0 },
	    {  1, -1,  0 } };

  
  
  // Lattice weights
  
  _omega = { 4.0 / 9.0,
	     1.0 / 9.0,
	     1.0 / 9.0,
	     1.0 / 9.0,
	     1.0 / 9.0,
	     1.0 / 36.0,
	     1.0 / 36.0,
	     1.0 / 36.0,
	     1.0 / 36.0 };


  // Reverse directions
  _reverse = { 0, 3, 4, 1, 2, 7, 8, 5, 6 };

  
    // // MRT matrix
    // _M.resize( 9 );
    // for(uint i = 0 ; i < 9 ; i++)
    // 	_M[i].resize(9);

    // _M[0] = { 1,  1,  1,  1,  1,  1,  1,  1,  1};
    // _M[1] = {-4, -1, -1, -1, -1,  2,  2,  2,  2};
    // _M[2] = { 4, -2, -2, -2, -2,  1,  1,  1,  1};
    // _M[3] = { 0,  1,  0, -1,  0,  1, -1, -1,  1};
    // _M[4] = { 0, -2,  0,  2,  0,  1, -1, -1,  1};
    // _M[5] = { 0,  0,  1,  0, -1,  1,  1, -1, -1};
    // _M[6] = { 0,  0, -2,  0,  2,  1,  1, -1, -1};
    // _M[7] = { 0,  1, -1,  1, -1,  0,  0,  0,  0};
    // _M[8] = { 0,  0,  0,  0,  0,  1, -1,  1, -1};


    // _invM.resize( 9 );
    // for(uint i = 0 ; i < 9 ; i++)
    // 	_invM[i].resize(9);

    // const double c1 = 1.0/9.0,
    // 	c2 = 1.0 / 6.0,
    // 	c3 = 1.0 / 4.0;

    // _invM[0] = {c1,   -c1,      c1,      0,     0,       0,     0,      0,    0};
    // _invM[1] = {c1,   -c1/4,   -c1/2,    c2,   -c2,      0,     0,      c3,   0};
    // _invM[2] = {c1,   -c1/4,   -c1/2,    0,     0,       c2,   -c2,    -c3,   0};
    // _invM[3] = {c1,   -c1/4,   -c1/2,   -c2,    c2,      0,     0,      c3,   0};
    // _invM[4] = {c1,   -c1/4,   -c1/2,    0,     0,      -c2,    c2,    -c3,   0};
    // _invM[5] = {c1,    c1/2,    c1/4,    c2,    c2/2,    c2,    c2/2,   0,    c3};
    // _invM[6] = {c1,    c1/2,    c1/4,   -c2,   -c2/2,    c2,    c2/2,   0,   -c3};
    // _invM[7] = {c1,    c1/2,    c1/4,   -c2,   -c2/2,   -c2,   -c2/2,   0,    c3};
    // _invM[8] = {c1,    c1/2,    c1/4,    c2,    c2/2,   -c2,   -c2/2,   0,   -c3};

    
}


// Default destructor
D2Q9::~D2Q9() {

  for(uint i = 0 ; i < _q ; i++)
    delete &_lvel[i];

  delete &_lvel;

  delete &_reverse;

}



/** Acces members */

const uint& D2Q9::d() const {
  
    return _d;
    
}
