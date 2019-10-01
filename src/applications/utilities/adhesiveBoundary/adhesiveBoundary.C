#include <iostream>

#include <iomanip>

#include <spotSampleCreator.H>

#include <dictionary.H>


using namespace std;



int main( int argc, char **argv ) {


    // Read boundaries from dictionary

    dictionary dict("properties/adhesiveProperties");

    vector<string> boundaries = dict.bracedEntriesNames("Boundaries");


    // Spot factory
    
    spotSampleCreator spcreator;
    
    
    // Move over boundaries and assign coefficients

    for( auto bdname : boundaries ) {

	string sptype = dict.lookUpOrDefault<string>("Boundaries/" + bdname + "/sampleType", "uniform");

	unique_ptr<spotSample> spots = spcreator.create(sptype, bdname);


	// Compute coefficients

	vector< pair<uint, scalar> > Gads = spots->computeAds();


	

	// Distribute over processors

	dictionary pdict("properties/parallel");

	uint np = (uint)pdict.lookUp<scalar>("numProc");

	
	ifstream infile;

	infile.open("lattice/points");

	uint npoints;

	infile >> npoints;

	infile.close();


	
	infile.open( ("lattice/lattice.graph.part." + to_string(np)).c_str()  );

	if( infile.is_open() == false ){
	
	    cout << " [ERROR]  Unable to find file lattice/lattice.graph.part." << np << endl;
	
	    exit(1);
	
	}	

	vector< map<uint, uint> > globalToLocal(np);

	vector<uint> pidIds(npoints);

	vector<uint> total(np);

	std::fill( total.begin(), total.end(), 0 );


	for( uint i = 0 ; i < npoints ; i++ ) {

	    uint pid;
	    
	    infile >> pid;

	    pidIds[i] = pid;

	    total[pid]++;
	    
	    globalToLocal[pid][i] = total[pid] - 1;

	}
       
	infile.close();
	

	ofstream outfile;
	
	for( uint pid = 0 ; pid < np ; pid++ ) {

	    outfile.open( ("processor" + to_string(pid) + "/lattice/" + bdname + "_gads").c_str()  );


	    for( auto g : Gads ) {

		if(pidIds[g.first] == pid) {

		    outfile << globalToLocal[pid][g.first] << " " << g.second << endl;

		}

	    }
	    
	    
	    outfile.close();

	}
	
    }


}
