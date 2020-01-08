#ifndef MLEVELKMETISDECOMP_H
#define MLEVELKMETISDECOMP_H

#include <vector>

#include "../meshInclude/basicMesh.H"


void mlevelKmetisDecomp( std::vector<uint>& owner, basicMesh& mesh, uint np );


#endif // MLEVELKMETISDECOMP_H
