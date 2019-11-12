import numpy as np

def check_active_points(vtkCells, npoints):
    """
    Check if point is considered according to cell arrays
    """

    active = np.zeros( (npoints,1), dtype=np.int64 )

    for cell in vtkCells:

        for c in cell:

            active[c] = 1
                               

    return active
    
