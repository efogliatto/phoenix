#include <timeOptions.H>

#include <scalarVector.H>


using namespace std;

int main() {

    dictionary dict("properties/macroProperties");

    cout << dict.lookUp<string>("EOS/model") << endl;

    cout << dict.lookUp<scalar>("EOS/a_vdw") << endl;

    cout << dict.lookUp<int>("EOS/G") << endl;


    scalarVector V( dict.lookUp< vector<scalar> >("f/Lambda") );

    cout << V << endl;



    scalarVector A( dict.lookUp< vector<scalar> >("g/Lambda") );

    cout << A << endl;



    timeOptions Time;


}
