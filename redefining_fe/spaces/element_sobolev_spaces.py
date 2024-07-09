from firedrake import *
import numpy as np
import sympy as sp
from ufl.sobolevspace import SobolevSpace
from redefining_fe.spaces.polynomial_spaces import PolynomialSpace
import matplotlib.pyplot as plt


class ElementSpaceTriple():

    def __init__(self, v, sobolev, interp_domain):
        assert isinstance(v, PolynomialSpace)
        assert isinstance(sobolev, SobolevSpace)
        # assert isinstance(interp_domain, ContinuousFunctionSpace)

        self.v = v
        self.w = sobolev
        self.w_I = interp_domain


class ElementSobolevSpace(SobolevSpace):
    """
    Representation of a Sobolev space on a single cell

    :param: *underlying_space*: The UFL representation of the Sobolev Space
    :param: *domain*: (Optional) the cell defined over- if not originally provided it should be provided during use.

    """

    def __init__(self, underlying_space, domain=None):
        self.domain = domain
        super(ElementSobolevSpace, self).__init__(underlying_space.name,
                                                  underlying_space.parents)

    def plot(self, ax, coord, trace_entity, g, **kwargs):
        raise NotImplementedError("Plotting not implemented for", str(type(self)))

    def pullback(self, v, trace_entity, g):
        raise NotImplementedError("Pullback not implemented for", str(type(self)))


class CellH1(ElementSobolevSpace):

    def __init__(self, cell):
        super(CellH1, self).__init__(H1, cell)

    def pullback(self, v, trace_entity, g):
        return v

    def plot(self, ax, coord, trace_entity, g, **kwargs):
        # plot dofs of the type associated with this space
        ax.scatter(*coord, **kwargs)

    def __repr__(self):
        return "H1"


class CellHDiv(ElementSobolevSpace):

    def __init__(self, cell):
        super(CellHDiv, self).__init__(HDiv, cell)

    def pullback(self, v, trace_entity, g):
        entityBasis = np.array(trace_entity.basis_vectors())
        cellEntityBasis = np.array(self.domain.basis_vectors(entity=trace_entity))
        basis = np.matmul(entityBasis, cellEntityBasis)

        def apply(*x):
            if len(v(*x)) == 2:
                result = np.cross(np.array(v(*x)).squeeze(), basis)
                # vec = np.matmul(np.array([[0, 1], [-1, 0]]), basis.T)
            elif trace_entity.dimension == 2:
                result = np.dot(np.array(v(*x)), np.cross(basis[0], basis[1]))
            else:
                raise ValueError("Immersion of HDiv edges not defined in 3D")
            if isinstance(result, np.float64):
                # todo: might always be a float
                return (result,)
            return tuple(result)
        return apply

    def plot(self, ax, coord, trace_entity, g, **kwargs):
        # plot dofs of the type associated with this space
        entityBasis = np.array(trace_entity.basis_vectors())
        cellEntityBasis = np.array(self.domain.basis_vectors(entity=trace_entity))
        basis = np.matmul(entityBasis, cellEntityBasis)
        if len(coord) == 2:
            vec = vec = np.matmul(np.array([[0, 1], [-1, 0]]), basis.T)
        else:
            vec = np.cross(basis[0], basis[1])
        ax.quiver(*coord, *vec, **kwargs)

    def __repr__(self):
        return "HDiv"


class CellHCurl(ElementSobolevSpace):

    def __init__(self, cell):
        super(CellHCurl, self).__init__(HCurl, cell)

    def pullback(self, v, trace_entity, g):
        tangent = np.array(trace_entity.basis_vectors())
        subEntityBasis = np.array(self.domain.basis_vectors(entity=trace_entity))

        def apply(*x):
            # breakpoint()
            result = np.dot(np.matmul(tangent, subEntityBasis),
                            np.array(v(*x)))
            if isinstance(result, np.float64):
                return (result,)
            return tuple(result)
        return apply

    def plot(self, ax, coord, trace_entity, g, **kwargs):
        # plot dofs of the type associated with this space
        tangent = np.array(trace_entity.basis_vectors())
        subEntityBasis = np.array(self.domain.basis_vectors(entity=trace_entity))
        vec = np.matmul(tangent, subEntityBasis)[0]
        ax.quiver(*coord, *vec, **kwargs)

    def __repr__(self):
        return "HCurl"


class CellH2(ElementSobolevSpace):

    def __init__(self, cell):
        super(CellH2, self).__init__(H2, cell)

    def pullback(self, v, trace_entity, g):
        # Compute grad v and then dot with tangent rotated according to the group member
        tangent = np.array(g(np.array(self.domain.basis_vectors())[0]))

        def apply(*x):
            X = sp.DeferredVector('x')
            dX = tuple([X[i] for i in range(self.domain.dim())])
            compute_v = v(*dX, sym=True)
            grad_v = [sp.diff(compute_v, dX[i]) for i in range(len(dX))]
            eval_grad_v = [comp.evalf(subs=dict(zip(dX, v.attach_func(*x)))) for comp in grad_v]
            result = np.dot(tangent, np.array(eval_grad_v))

            if not hasattr(result, "__iter__"):
                return (result,)
            return tuple(result)
        return apply

    def plot(self, ax, coord, trace_entity, g, **kwargs):
        circle1 = plt.Circle(coord, 0.075, fill=False, **kwargs)
        ax.add_patch(circle1)

    def __repr__(self):
        return "H2"


class CellH3(ElementSobolevSpace):

    def __init__(self, cell):
        super(CellH3, self).__init__(H2, cell)

    def pullback(self, v, trace_entity, g):
        b0, b1 = self.domain.basis_vectors()
        tangent0 = np.array(g(b0))
        tangent1 = np.array(g(b1))

        def apply(*x):
            X = sp.DeferredVector('x')

            dX = tuple([X[i] for i in range(self.domain.dim())])
            hess_v = [[sp.diff(v(*dX, sym=True), dX[i], dX[j]) for i in range(len(dX))] for j in range(len(dX))]
            print(hess_v)
            eval_grad_v = [[hess_v[i][j].evalf(subs=dict(zip(dX, v.attach_func(*x)))) for i in range(len(dX))] for j in range(len(dX))]
            result = np.dot(np.matmul(tangent0, np.array(eval_grad_v)), tangent1)
            if not hasattr(result, "__iter__"):
                return (result,)
            return tuple(result)
        return apply

    def plot(self, ax, coord, trace_entity, g, **kwargs):
        circle1 = plt.Circle(coord, 0.15, fill=False, **kwargs)
        ax.add_patch(circle1)

    def __repr__(self):
        return "H3"


class CellL2(ElementSobolevSpace):

    def __init__(self, cell):
        super(CellL2, self).__init__(L2, cell)

    def pullback(self, v, trace_entity):
        return 0

    def plot(self, ax, coord, trace_entity, g, **kwargs):
        # plot dofs of the type associated with this space
        ax.scatter(*coord, **kwargs)

    def __repr__(self):
        return "L2"
