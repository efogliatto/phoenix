#include <lbEquation.H>

/** Default constructor */

lbEquation::lbEquation( const std::string& name ) : ename (name) {}


/** Default destructor */

lbEquation::~lbEquation() {}



/** Collision process */

const void lbEquation::collision() {}


/** Streamming process */

const void lbEquation::streamming() {}


/** Update macroscopic field */

const void lbEquation::updateMacroField() {}
