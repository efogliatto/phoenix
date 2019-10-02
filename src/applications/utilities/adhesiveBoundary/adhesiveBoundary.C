#include <iostream>

#include <iomanip>

#include <spotSampleCreator.H>

#include <dictionary.H>

#include <localIndexer.H>


using namespace std;



int main( int argc, char **argv ) {


    // Read boundaries from dictionary

    dictionary dict("properties/adhesiveProperties");

    vector<string> boundaries = dict.bracedEntriesNames("Boundaries");


    // Spot factory
    
    spotSampleCreator spcreator;
    
    
    // Move over boundaries and assign coefficients

    for( auto bdname : boundaries ) {


	// Create spots

	string sptype = dict.lookUpOrDefault<string>("Boundaries/" + bdname + "/sampleType", "uniform");

	unique_ptr<spotSample> spots = spcreator.create(sptype, bdname);


	// Compute coefficients

	vector< pair<uint, scalar> > Gads = spots->computeAds();

	

	// Distribute over processors

	localIndexer indexer;

	
	for( uint pid = 0 ; pid < indexer.np() ; pid++ ) {


	    // First check total number of boundary nodes for this processor

	    uint count(0);

	    for( const auto& g : Gads ) {

		int lid = indexer.globalToLocal( g.first, pid );

		if(lid != -1)
		    count++;

	    }


	    if( count > 0 ) {
	    
		ofstream outfile;
	    
		outfile.open( ("processor" + to_string(pid) + "/lattice/" + bdname + "_gads").c_str()  );

		outfile << count << endl;

		for( const auto& g : Gads ) {

		    int lid = indexer.globalToLocal( g.first, pid );

		    if( lid != -1 )
			outfile << lid << "  " << g.second << endl;

		}


	    }


	}
	
	
    }


}
