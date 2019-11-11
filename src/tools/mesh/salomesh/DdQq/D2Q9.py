import numpy as np


class D2Q9:

    """
    D2Q9 lattice model
    """

    def __init__(self):

        """
        Constructor
        """

        # DdQq constructor

        super().__init__()
        

        # Velocity set
        
        self.__vel = np.array( [[0,0,0], [1,0,0], [0,1,0], [-1,0,0], [0,-1,0], [1,1,0], [-1,1,0], [-1,-1,0], [1,-1,0] ], dtype=np.int64)


        # Reverse indices

        self.__rev = np.array( [0, 3, 4, 1, 2, 7, 8, 5, 6], dtype=np.int64 )


        pass


    # Return velocities

    def velocities(self):
        """
        Velocity set
        """
        return self.__vel


    # Return reverse indices

    def reverse(self):
        """
        Reverse indices
        """
        return self.__rev    
