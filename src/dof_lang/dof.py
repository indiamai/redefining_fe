
from cell_complex.cells import Point, Edge


class Pairing():

    def __init__(self, entity, space):
        self.entity = entity
        self.space = space


class DeltaPairing(Pairing):

    def __init__(self, entity, space):
        super(DeltaPairing, self).__init__(entity, space)

    def __call__(self, kernel, v):
        assert isinstance(kernel, DeltaKernel)
        return v(*kernel.pt)

    def __repr__(self, kernel):
        return "v(%s)" % str(kernel)
    
    def immerse(self, E):
        return DeltaPairing(E, self.space)


class L2InnerProd(Pairing):
    """ need to think about the abstraction level here - 
    are we wanting to define them as quadrature now? or defer this?
    """
    def __init__(self, entity, space):
        super(L2InnerProd, self).__init__(entity, space)

    def __call__(self, x, v):
        # evaluates integral
        pass

    def __repr__(self, kernel):
        return "integral_{1}({0} * v dx)".format(str(kernel), str(self.entity))
    
    def immerse(self, E):
        raise NotImplementedError("how does the cell transfor")


class DeltaKernel():

    def __init__(self, x):
        self.pt = x

    def __repr__(self):
        x = list(map(str, list(self.pt)))
        return ','.join(x)
    
    def immerse(self, attachment, E):
        return DeltaKernel(attachment(self.pt))


class TangentKernel():

    def __init__(self, E):
        self.tangent = E.basis_vectors(return_coords=True)[0]

    def __repr__(self):
        return str(self.tangent)
    
    def immerse(self, attachment, E):
        return TangentKernel(E)


class DOF():

    def __init__(self, pairing, kernel):
        self.pairing = pairing
        self.kernel = kernel

    def __call__(self, fn):
        return self.pairing(self.kernel, fn)

    def __repr__(self):
        return self.pairing.__repr__(self.kernel)
    
    def trace(self, attachment, cell):
        return DOF(self.pairing.immerse(cell),
                   self.kernel.immerse(attachment, cell))


def construct_point_eval(x, E, V):
    delta_x = DeltaKernel(x)
    return DOF(DeltaPairing(E, V), delta_x)


def construct_tangent_dof(E, V):
    return DOF(L2InnerProd(E, V), TangentKernel(E))
