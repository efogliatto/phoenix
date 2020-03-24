#ifndef LINEARSCALARFIELD_H
#define LINEARSCALARFIELD_H

#include <scalarField.H>

#include <latticeMesh.H>

void linearScalarField( scalarField& field, const latticeMesh& mesh, const std::string& fname, const std::string& sptype );

#endif // LINEARSCALARFIELD_H
