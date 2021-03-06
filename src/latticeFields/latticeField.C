#include <latticeField.H>

/** Read field using ensight format */

const void latticeField::read() {}


/** Default constructor */

latticeField::latticeField( const latticeMesh& m, timeOptions& t, const std::string& nm ) : name(nm), mesh(m), Time(t)  {}


/** Default destructor */

latticeField::~latticeField() {}



// Acces members. Forced interface

/** Synchronization across procceses */

const void latticeField::sync() {}


/** Start sync */

const void latticeField::startSync() {}


/** End sync */

const void latticeField::endSync() {}



/** Write field using ensight format */

const void latticeField::write() const {}
