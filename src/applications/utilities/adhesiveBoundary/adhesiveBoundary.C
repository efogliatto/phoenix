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

	for(auto g : Gads)
	    cout << g.first << "  " << g.second << endl;

    }


}
