import numpy as np


def vtk_cells(grid, lmodel):


    
    nx = grid[0]
    ny = grid[1]
    nz = grid[2]
    

    
    if lmodel.D() == 2:

        vtkCells = np.zeros(  ( (nx-1)*(ny-1),4)  )

        ncells = 0

    
        for j in range(ny-1):

            for i in range(nx-1):    

    	        vtkCells[ncells,0] = i + j*nx
    	        vtkCells[ncells,1] = i + j*nx + 1
    	        vtkCells[ncells,2] = i + j*nx + nx
    	        vtkCells[ncells,3] = i + j*nx + nx + 1
	
    	        ncells = ncells + 1


                


    elif lmodel.D() == 3:

        vtkCells = np.zeros(  ( (nx-1)*(ny-1)*(nz-1),8)  )    

        ncells = 0

        for k in range(nz-1):
    
            for j in range(ny-1):
	
                for i in range(nx-1):

                    vtkCells[ncells,0] = i + j*nx + k*nx*ny

                    vtkCells[ncells,1] = i  +  j*nx + 1  +  k*nx*ny

                    vtkCells[ncells,2] = i  +  j*nx + nx  +  k*nx*ny

                    vtkCells[ncells,3] = i  +  j*nx + nx + 1  +  k*nx*ny

                    vtkCells[ncells,4] = i  +  j*nx  +  (k+1)*nx*ny

                    vtkCells[ncells,5] = i  +  j*nx + 1  +   (k+1)*nx*ny

                    vtkCells[ncells,6] = i  +  j*nx + nx  +  (k+1)*nx*ny

                    vtkCells[ncells,7] = i  +  j*nx + nx + 1  +  (k+1)*nx*ny

                    ncells = ncells + 1
    

                    
            
    return vtkCells
