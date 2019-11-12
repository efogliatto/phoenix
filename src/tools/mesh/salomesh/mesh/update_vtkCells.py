import numpy as np

def update_vtkCells(vtkCells, oldToNew):
    """
    Update vtk cell using new point indices
    """

    sp = vtkCells.shape
    
    newVtk = np.zeros( sp, dtype=np.int64 )

    for i in range(sp[0]):

        for j in range(sp[1]):

            newVtk[i,j] = oldToNew[ vtkCells[i,j] ]


    return newVtk
