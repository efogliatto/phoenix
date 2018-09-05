#ifndef BASICMESH_H
#define BASICMESH_H

#include <dataTypes.H>

#include <latticeModelCreator.H>

#include <vector>



class basicMesh {

protected:

    /** Total number of points per patch */

    uint nPoints;


    /** Points coordinates */

    std::vector< std::vector<int> > points;


    /** Neighbour indices array */

    int** nb;


    /** Lattice model */

    const latticeModel* lbm;



public:


    /** Default constructor */

    basicMesh( const latticeModel* model );


    /** Default destructor */

    ~basicMesh();


};

#endif // BASICMESH_H
