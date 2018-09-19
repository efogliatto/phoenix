#include <timeOptions.H>

#include <scalarVector.H>


using namespace std;

int main() {

    dictionary dict("properties/macroProperties");

    vector<string> entry = dict.bracedEntriesNames("Navier-Stokes");

    for(auto s : entry)
	cout << s << endl;


    return 0;
    
}
