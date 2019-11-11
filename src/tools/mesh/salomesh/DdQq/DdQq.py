from .D2Q9 import D2Q9
from .D3Q15 import D3Q15

def DdQq( model="D2Q9" ):

    """
    Model set
    """

    if model == "D2Q9":
        return D2Q9()

    elif model == "D3Q15":
        return D3Q15()

    else:
        assert False, "Unrecognized lattice model"
