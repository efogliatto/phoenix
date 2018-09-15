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


    /** Number of local points */

    const uint local() const { return parallel.local(); }       
    
    
    /** Nodes per patch */

    const uint npp( const uint& id ) const { return parallel.npp(id); }

    
    /** Reference to model */

    const latticeModel* lmodel() const { return lbmodel; }


    /** Reference to neighbours array */

    const std::vector< std::vector<int> >& nbArray() const { return nb; }


    /** Reference to boundary indiced */

    const std::vector<uint>& boundaryNodes( const std::string& bdname ) const;

    const std::map< std::string, std::vector<uint> >& boundaries() const { return boundary; }
    
    
};

#endif // LATTICEMESH_H
