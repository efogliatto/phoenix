#ifndef KMETISDECOMP_H
#define KMETISDECOMP_H

#include <vector>

#include "../meshInclude/basicMesh.H"


void kmetisDecomp( std::vector<uint>& owner, basicMesh& mesh, uint np, uint maxneigh );


#endif // KMETISDECOMP_H
