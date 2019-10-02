#include <localIndexer.H>

/** Default constructor */

localIndexer::localIndexer() {}


/** Constructor with vector of indices */

localIndexer::localIndexer( const std::vector<uint>& pid ) {}


/** Destructor */

localIndexer::~localIndexer() {}


/** Global to local index */

const std::pair<int, int> localIndexer::globalToLocal( const uint id ) const {}
