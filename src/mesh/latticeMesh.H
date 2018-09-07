#ifndef LATTICEMESH_H
#define LATTICEMESH_H

#include <vector>

#include <mpiInfo.H>


class latticeMesh {

    
protected:


    /** Total number of points per patch */

    uint nPoints;


    /** Points coordinates */

    std::vector< std::vector<int> > points;


    /** Neighbour indices array */

    int** nb;
    

    /** MPI comm information */

    const mpiInfo parallel;


    /** Read mesh points */

    const void readPoints();



public:

    
    /** Default constructor */

    latticeMesh( const int& pid );


    /** Default destructor */

    ~latticeMesh();

};

#endif // LATTICEMESH_H
