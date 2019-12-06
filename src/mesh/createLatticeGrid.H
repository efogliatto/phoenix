#ifndef CREATELATTICEGRID_H
#define CREATELATTICEGRID_H

#include <dataTypes.H>

#include <vector>

#include <latticeModelCreator.H>


void createLatticeGrid( std::vector< std::vector<uint> >& points, const latticeModel* lbmodel, const uint nx, const uint ny, const uint nz = 0  );


#endif // CREATELATTICEGRID_H
