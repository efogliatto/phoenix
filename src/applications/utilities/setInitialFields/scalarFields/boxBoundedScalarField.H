#ifndef BOXBOUNDEDSCALARFIELD_H
#define BOXBOUNDEDSCALARFIELD_H

#include <scalarField.H>

#include <latticeMesh.H>

void boxBoundedScalarField( scalarField& field, const latticeMesh& mesh, const std::string& fname, const std::string& sptype );

#endif // BOXBOUNDEDSCALARFIELD_H
