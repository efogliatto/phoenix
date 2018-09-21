#include <heatSource.H>

using namespace std;



/** Constructor */

heatSource::heatSource( const string& dictName,
			const string& eqName,
			const latticeMesh& mesh,
			timeOptions& Time )

    : _mesh(mesh),
      _source(mesh, Time, "HeatSource", IO::NO_READ, IO::NO_WRITE) {    

}



/** Destructor */

heatSource::~heatSource() {}
