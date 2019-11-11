import numpy as np


class D3Q15:

    """
    D3Q15 lattice model
    """

    def __init__(self):

        """
        Constructor
        """

        # DdQq constructor

        super().__init__()
        

        # Velocity set
        
        self.__vel = np.array( [ [  0,  0,  0 ],
                                 [  1,  0,  0 ],
                                 [ -1,  0,  0 ],
                                 [  0,  1,  0 ],
                                 [  0, -1,  0 ],
                                 [  0,  0,  1 ],
                                 [  0,  0, -1 ],
                                 [  1,  1,  1 ],
                                 [ -1,  1,  1 ],
                                 [  1, -1,  1 ],
                                 [ -1, -1,  1 ],
                                 [  1,  1, -1 ],
                                 [ -1,  1, -1 ],
                                 [  1, -1, -1 ],
                                 [ -1, -1, -1 ] ], dtype=np.int64)

        # Reverse indices

        self.__rev = np.array( [0, 2, 1, 4, 3, 6, 5, 14, 13, 12, 11, 10, 9, 8, 7], dtype=np.int64 )
        
        
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


    # Return name

    def name(self):
        """
        Model name
        """
        return "D3Q15"    


    # Dimension

    def D(self):
        """
        Lattice dimension
        """
        return np.int64(3)    
