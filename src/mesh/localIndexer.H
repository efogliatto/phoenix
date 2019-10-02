#ifndef LOCALINDEXER_H
#define LOCALINDEXER_H

#include <dictionary.H>

class localIndexer {


protected:

    /** Global to local list */

    std::vector<  std::vector< std::pair<int, int> >  > _idx;
    
    
public:

    /** Default constructor */

    localIndexer();


    /** Constructor with vector of indices */

    localIndexer( const std::vector<uint>& pid );


    /** Destructor */

    ~localIndexer();


    /** Global to local index */

    const std::pair<int, int> globalToLocal( const uint id ) const;

};

#endif // LOCALINDEXER_H
