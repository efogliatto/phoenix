#ifndef BGSPHEREBOUNDEDSCALARFIELD_H
#define BGSPHEREBOUNDEDSCALARFIELD_H

#include <scalarField.H>

#include <latticeMesh.H>

void bgSphereBoundedScalarField( scalarField& field, const latticeMesh& mesh, const std::string& fname, const std::string& sptype );

#endif // BGSPHEREBOUNDEDSCALARFIELD_H
