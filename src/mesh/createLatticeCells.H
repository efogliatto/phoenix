#ifndef CREATELATTICECELLS_H
#define CREATELATTICECELLS_H

#include <dataTypes.H>

#include <vector>

#include <latticeModelCreator.H>


void createLatticeCells( std::vector< std::vector<uint> >& cells, const latticeModel* lbmodel, const uint nx, const uint ny, const uint nz = 0  );

#endif // CREATELATTICECELLS_H
