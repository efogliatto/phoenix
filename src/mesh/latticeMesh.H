#ifndef LATTICEMESH_H
#define LATTICEMESH_H

#include <vector>

#include <mpiInfo.H>

#include <latticeModelCreator.H>

#include <map>


class latticeMesh {

    
protected:


    /** Total number of points per patch */

    uint nPoints;


    /** Points coordinates */

    std::vector< std::vector<int> > points;


    /** Neighbour indices array */

    std::vector< std::vector<int> > nb;
    

    /** MPI comm information */

    const mpiInfo parallel;

    
    /** Lattice Model */

    latticeModel* lbmodel;


    /** Boundary nodes */

    std::map< std::string, std::vector<uint> > boundary;



    

    


    /** Read mesh points */

    const void readPoints();


    /** Read neighbours */

    const void readNeighbours();

    
    /** Read boundary nodes */

    const void readBoundaryNodes();
    



public:

    
    /** Default constructor */

    latticeMesh( const int& pid, const bool& msg = true );


    /** Default destructor */

    ~latticeMesh();


    /** Process id */

    const int pid() const { return parallel.id(); }

    
    /** Total number of processors */

    const int wsize() const { return parallel.wsize(); }
    

    /** Total number of points */

    const uint npoints() const { return nPoints; }

    
    /** Nodes per patch */

    const uint npp( const uint& id ) const { return parallel.npp(id); }
    
    
};

#endif // LATTICEMESH_H
