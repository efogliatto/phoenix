#include <iostream>
#include <algebra.H>


using namespace std;

int main() {

    const scalarMatrix M( { { 1,  1,  1,  1,  1,  1,  1,  1,  1},
	                    {-4, -1, -1, -1, -1,  2,  2,  2,  2},
			    { 4, -2, -2, -2, -2,  1,  1,  1,  1},
			    { 0,  1,  0, -1,  0,  1, -1, -1,  1},
			    { 0, -2,  0,  2,  0,  1, -1, -1,  1},
			    { 0,  0,  1,  0, -1,  1,  1, -1, -1},
			    { 0,  0, -2,  0,  2,  1,  1, -1, -1},
			    { 0,  1, -1,  1, -1,  0,  0,  0,  0},
			    { 0,  0,  0,  0,  0,  1, -1,  1, -1}  } );


    const scalarVector V( M[1] );

    scalarVector res( M[1] );
    
    cout << M << endl;
    
    cout << V << endl;

    M.matDotVec(V, res);

    cout << res << endl << endl;



    vector<scalar> c = {0,0,0};

    sparseScalarMatrix S( {1,2,3} );

    S.addElement(4,1,2);

    S.matDotVec( {1,1,1} ,c);

    for (auto i : c)
	cout << i << endl;

    


}
