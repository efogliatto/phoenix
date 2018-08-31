#include <latticeModel.H>

using namespace std;



/* ----------------------  Public member functions ----------------------  */



// Constructors and destructors

/** Default constructor */

latticeModel::latticeModel() {}


/** Default destructor */

latticeModel::~latticeModel() {}





// Access members

/** Discrete lattice velocities */

const vector< vector<int> >& latticeModel::lvel() const {

  return _lvel;

}


/** Dimension */

const uint& latticeModel::d() const {

  return _d;
  
}


/** Number of lattice velocities */

const uint latticeModel::q() const {
  
  return _q;
    
}


/** Square of sound speed (normaliced by lattice velocity) */

const scalar& latticeModel::cs2() const {

  return _cs2;
  
}


/** Lattice weights */

const vector<scalar>& latticeModel::omega() const {
  
  return _omega;
    
}


/** Reverse indices */

const vector<uint>& latticeModel::reverse() const {
  
  return _reverse;
  
}
