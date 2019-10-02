#ifndef LOCALINDEXER_H
#define LOCALINDEXER_H

#include <dictionary.H>

#include <map>

class localIndexer {


protected:

    /** Global to local list : (pid, node) */

    std::vector<  std::map< uint, uint>  > _idx;

    
    /** Local to global list : (local, global) */

    std::vector<  std::map< uint, uint>  > _l2g;    

    
    /** Total number of nodes per processor */

    std::vector<uint> nodesPerProc;
    
    
    /** Construct global to local list */

    const void createIdxList( const std::vector<uint>& pids, const std::vector< std::vector<int> >& nb );
   
    
    
public:

    /** Default constructor */

    localIndexer();


    /** Constructor with vector of indices */

    localIndexer( const std::vector<uint>& pids, const std::vector< std::vector<int> >& nb );


    /** Destructor */

    ~localIndexer();


    /** Global to local index */

    const int globalToLocal( const uint id, const uint pid ) const;


    /** Local to global index */

    const int localToGlobal( const uint id, const uint pid ) const;


    /** Local to local index */

    const int localToLocal( const uint id, const uint ownPid, const uint otherPid ) const;


    /** Total number of processors */

    const uint np() const;

};

#endif // LOCALINDEXER_H
