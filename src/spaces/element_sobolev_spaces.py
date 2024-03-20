from firedrake import *
import numpy as np
from ufl.sobolevspace import SobolevSpace
from spaces.polynomial_spaces import PolynomialSpace

def apply_pullback(*x):
    result = np.dot(np.matmul(tangent, subEntityBasis), np.array(v(*x)))
    if isinstance(result, np.float64):
        return (result,)
    return tuple(result)

class ElementSpaceTriple():

    def __init__(self, v, sobolev, interp_domain):
        assert isinstance(v, PolynomialSpace)
        assert isinstance(sobolev, SobolevSpace)
        # assert isinstance(interp_domain, ContinuousFunctionSpace)

        self.v = v
        self.w = sobolev
        self.w_I = interp_domain


class ElementSobolevSpace(SobolevSpace):

    def __init__(self, underlying_space, domain=None):
        self.domain = domain
        super(ElementSobolevSpace, self).__init__(underlying_space.name,
                                                  underlying_space.parents)

    def trace(self, v):
        raise NotImplementedError("Trace not implemented for", str(type(self)))

    def pullback(self, v, trace_entity):
        raise NotImplementedError("Pullback not implemented for", str(type(self)))


class CellH1(ElementSobolevSpace):

    def __init__(self, cell):
        super(CellH1, self).__init__(H1, cell)

    def pullback(self, v, trace_entity):
        # temporarily everything is reference space?
        return v
    
    def __repr__(self):
        return "H1"


class CellHDiv(ElementSobolevSpace):

    def __init__(self, cell):
        super(CellHDiv, self).__init__(HDiv, cell)

    def pullback(self, v, trace_entity):
        entityBasis = np.array(trace_entity.basis_vectors())
        cellEntityBasis = np.array(self.domain.basis_vectors(entity=trace_entity))
        
        def apply(*x):
            result = np.cross(np.array(v(*x)),
                              np.matmul(entityBasis, cellEntityBasis))
            if isinstance(result, np.float64):
                return (result,)
            return tuple(result)
        return apply

    def __repr__(self):
        return "HDiv"


class CellHCurl(ElementSobolevSpace):

    def __init__(self, cell):
        super(CellHCurl, self).__init__(HCurl, cell)

    def pullback(self, v, trace_entity):
        tangent = np.array(trace_entity.basis_vectors())
        subEntityBasis = np.array(self.domain.basis_vectors(entity=trace_entity))

        def apply(*x):
            result = np.dot(np.matmul(tangent, subEntityBasis),
                            np.array(v(*x)))
            if isinstance(result, np.float64):
                return (result,)
            return tuple(result)
        return apply

    def __repr__(self):
        return "HCurl"


class CellL2(ElementSobolevSpace):

    def __init__(self, cell):
        super(CellL2, self).__init__(L2, cell)

    def pullback(self, v, trace_entity):
        return 0
    
    def __repr__(self):
        return "L2"
