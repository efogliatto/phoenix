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
    
    cout << M << endl;
    
    cout << V << endl;

    cout << M * V ;


}
