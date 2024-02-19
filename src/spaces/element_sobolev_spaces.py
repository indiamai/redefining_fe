from firedrake import *
from ufl.sobolevspace import SobolevSpace
from polynomial_spaces import PolynomialSpace


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

    def pullback(self, v):
        raise NotImplementedError("Pullback not implemented for", str(type(self)))


class CellH1(ElementSobolevSpace):

    def __init__(self, cell):
        super(CellH1, self).__init__(H1, cell)

    def pullback(self, v):
        # temporarily everything is reference space?
        return v


class CellHDiv(ElementSobolevSpace):

    def __init__(self, cell):
        super(CellHDiv, self).__init__(HDiv, cell)

    def pullback(self, v):
        # temporarily everything is reference space?
        return v


class CellL2(ElementSobolevSpace):

    def __init__(self, cell):
        super(CellL2, self).__init__(L2, cell)

    def pullback(self, v):
        # temporarily everything is reference space?
        return v
