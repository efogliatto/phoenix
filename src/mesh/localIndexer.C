#include <localIndexer.H>

#include <algorithm>


using namespace std;


/** Default constructor */

localIndexer::localIndexer() {

       
    // Read elements from file

    ifstream infile;

    infile.open("lattice/points");

    uint npoints;

    infile >> npoints;

    infile.close();



    vector<uint> pids(npoints);
    
    vector< vector<int> > nb;



    // Distribute over processors

    dictionary pdict("properties/parallel");

    uint np = (uint)pdict.lookUp<scalar>("numProc");


    if(np == 1) {
	
	std::fill(pids.begin(), pids.end(), 0);

    }


    else {

	infile.open( ("lattice/lattice.graph.part." + to_string(np)).c_str()  );

	if( infile.is_open() == false ){
	
	    cout << " [ERROR]  Unable to find file lattice/lattice.graph.part." << np << endl;
	
	    exit(1);
	
	}
	
	for(uint i = 0 ; i < npoints ; i++) {

	    uint p;

	    infile >> p;

	    pids[i] = p;

	}


    }

    
    createIdxList(pids, nb);



}


/** Constructor with vector of indices */

localIndexer::localIndexer( const vector<uint>& pids, const vector< vector<int> >& nb ) {

    createIdxList(pids, nb);

}


/** Destructor */

localIndexer::~localIndexer() {}


/** Global to local index */

const int localIndexer::globalToLocal( const uint id, const uint pid ) const {

    int p(-1);

    if( _idx[id].find(pid) != _idx[id].end() )
    	p = _idx[id].at(pid);
    
    return p;

}


/** Local to global index */

const int localIndexer::localToGlobal( const uint id, const uint pid ) const {

    int p(-1);
    
    
    if( pid < _l2g.size() ) {

	if( _l2g[pid].find(id) != _l2g[pid].end() )
	    p = _l2g[pid].at(id);

    }

    else {

	cout << " [ERROR]  Incorrect number of processors" << endl;
	
	exit(1);

    }


    return p;

}


/** Local to local index */

const int localIndexer::localToLocal( const uint id, const uint ownPid, const uint otherPid ) const {

    return -1;

}



/** Construct global to local list */

const void localIndexer::createIdxList( const vector<uint>& pids, const vector< vector<int> >& nb ) {


    // Total number of processors

    nodesPerProc.resize( *std::max_element(pids.begin(), pids.end()) + 1 );

    std::fill( nodesPerProc.begin(), nodesPerProc.end(), 0 );
    

    // Resize local indices

    _idx.resize( pids.size() );
    
    
    // Move over list and assign

    for( uint i = 0 ; i < pids.size() ; i++ ) {

    	uint p = pids[i];
	
    	_idx[i][p] = nodesPerProc[p];

    	nodesPerProc[p]++;	

    }



    // Create local to global list

    _l2g.resize( nodesPerProc.size() );

    for( uint i = 0 ; i < _idx.size() ; i++ ) {

    	for( const auto &p : _idx[i] ) {

    	    _l2g[p.first][p.second] = i;

    	}

    }
    

}


/** Total number of processors */

const uint localIndexer::np() const {

    return nodesPerProc.size();

}
