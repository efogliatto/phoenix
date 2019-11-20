#include <dictionary.H>

#include "readVTKInfo.H"

vtkInfo readVTKInfo() {


    vtkInfo vtk;

    uint status = 0;

    
    status = lookUpStringList( "start/initialFields", "scalarFields", &vtk.scalarFields, &vtk.nscalar);

    if( !status ) { vtk.nscalar = 0; }

    
    status = lookUpStringList( "start/initialFields", "vectorFields", &vtk.vectorFields, &vtk.nvector);

    if( !status ) { vtk.nvector = 0; }
    

    status = lookUpStringList( "start/initialFields", "pdfFields", &vtk.pdfFields, &vtk.npdf);

    if( !status ) { vtk.npdf = 0; }


    
    return vtk;

}
