#ifndef MPMETISDECOMP_H
#define MPMETISDECOMP_H

#include <vector>

#include "../meshInclude/basicMesh.H"


void mpmetisDecomp( std::vector<uint>& owner, basicMesh& mesh, uint np );


#endif // MPMETISDECOMP_H
