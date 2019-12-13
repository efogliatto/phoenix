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


/** Index matching velocity */

const int latticeModel::velocityIndex( const int i, const int j, const int k ) const {


    int idx(-1);


    for( uint vid = 0 ; vid < _q ; vid++ ) {

	if( _lvel[vid][0] == i ) {

	    if( _lvel[vid][1] == j ) {

	    	if( _lvel[vid][2] == k ) {

		    idx = vid;

		    vid = _q;

		}

	    }	    

	}

    }
    

    return idx;

}




/** Symmetric indices */

const vector<uint> latticeModel::symIdx( const string& symPlane ) const {

    vector<uint> idx;
    
    if( _symIdx.find( symPlane )  !=  _symIdx.end()  ) {
	
	idx = _symIdx.at( symPlane );

    }

    else {

	cout << endl << " [ERROR]  Unable to match symmetry plane " << symPlane << " with lattice model " << name() << endl;

	exit(1);

    }

    return idx;

}
