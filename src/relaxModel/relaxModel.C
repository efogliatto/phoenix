#include <relaxModel.H>

using namespace std;


/** Constructor */

relaxModel::relaxModel( const string& entry ) : _name("None") {}


/** Destructor */

relaxModel::~relaxModel() {}


/** Density-dependent relaxation factor */

const scalar relaxModel::tau( const scalar rho, const uint i ) { return 1; }
