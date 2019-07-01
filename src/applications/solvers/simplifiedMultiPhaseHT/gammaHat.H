#include <scalarField.H>

#include <vectorField.H>

#include <EOS.H>


const scalar gammaHat(const latticeMesh& mesh,
		      const scalarField& rho,
		      const scalarField& T,
		      const vectorField& U,
		      const scalar a1,
		      const scalar a2,
		      const scalar tau,
		      const EOS* eos,
		      const scalar Cv,
		      const uint id);
