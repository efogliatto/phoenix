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




    // MRT matrix
    
    M.setFromArray( { { 1,  1,  1,  1,  1,  1,  1,  1,  1},
	              {-4, -1, -1, -1, -1,  2,  2,  2,  2},
		      { 4, -2, -2, -2, -2,  1,  1,  1,  1},
		      { 0,  1,  0, -1,  0,  1, -1, -1,  1},
		      { 0, -2,  0,  2,  0,  1, -1, -1,  1},
		      { 0,  0,  1,  0, -1,  1,  1, -1, -1},
		      { 0,  0, -2,  0,  2,  1,  1, -1, -1},
		      { 0,  1, -1,  1, -1,  0,  0,  0,  0},
		      { 0,  0,  0,  0,  0,  1, -1,  1, -1} } );



    const scalar c1 = 1.0/9.0,
    	c2 = 1.0 / 6.0,
    	c3 = 1.0 / 4.0;
    
    
    invM.setFromArray( { {c1,   -c1,      c1,      0,     0,       0,     0,      0,    0},
	                 {c1,   -c1/4,   -c1/2,    c2,   -c2,      0,     0,      c3,   0},
			 {c1,   -c1/4,   -c1/2,    0,     0,       c2,   -c2,    -c3,   0},
			 {c1,   -c1/4,   -c1/2,   -c2,    c2,      0,     0,      c3,   0},
			 {c1,   -c1/4,   -c1/2,    0,     0,      -c2,    c2,    -c3,   0},
			 {c1,    c1/2,    c1/4,    c2,    c2/2,    c2,    c2/2,   0,    c3},
			 {c1,    c1/2,    c1/4,   -c2,   -c2/2,    c2,    c2/2,   0,   -c3},
			 {c1,    c1/2,    c1/4,   -c2,   -c2/2,   -c2,   -c2/2,   0,    c3},
			 {c1,    c1/2,    c1/4,    c2,    c2/2,   -c2,   -c2/2,   0,   -c3} } );


    // Symetric indices

    _symIdx["OX"] = {0,3,2,1,4,6,5,8,7};
    _symIdx["OY"] = {0,1,4,3,2,8,7,6,5};
    
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
