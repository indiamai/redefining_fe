from FIAT.polynomial_set import ONPolynomialSet
from FIAT.quadrature_schemes import create_quadrature
from FIAT import expansions, polynomial_set, reference_element
from itertools import chain
from redefining_fe.cells import CellComplexToFiat
from redefining_fe.utils import tabulate_sympy
import sympy as sp


class PolynomialSpace(object):
    """
    subdegree: the degree of the maximum degree Lagrange space that is spanned by this element. If this
    element's polynomial space does not include the constant function, this function should
    return -1.

    super degree: the degree of the minimum degree Lagrange space that spans this element.If this
    element contains basis functions that are not in any Lagrange space, this property should
    be None.

    Note that on a simplex cells, the polynomial space of Lagrange space is a complete polynomial
    space, but on other cells this is not true. For example, on quadrilateral cells, the degree 1
    Lagrange space includes the degree 2 polynomial xy.
    """

    def __init__(self, subdegree, superdegree):
        self.subdegree = subdegree
        self.superdegree = superdegree

    def complete(self):
        return self.subdegree == self.superdegree

    def to_ON_polynomial_set(self, ref_el):
        # how does super/sub degrees work here
        if not isinstance(ref_el, reference_element.Cell):
            ref_el = CellComplexToFiat(ref_el)
        return ONPolynomialSet(ref_el, self.subdegree)

    def __repr__(self):
        if self.complete():
            return "P" + str(self.subdegree)
        else:
            return "Psub" + str(self.subdegree) + "sup" + str(self.superdegree)

    def __mul__(self, x):
        if isinstance(x, sp.Symbol):
            return ConstructedPolynomialSpace([x], [self])
        else:
            raise TypeError(f'Cannot multiply a PolySpace with {type(x)}')

    __rmul__ = __mul__

    def __add__(self, x):
        return ConstructedPolynomialSpace([1, 1], [self, x])

    def restrict(self, min_degree, max_degree):
        return RestrictedPolynomialSpace(min_degree, max_degree)


class RestrictedPolynomialSpace(PolynomialSpace):
    """
    Represents a polynomial space of all polynomials between two degrees.

    :param: min_degree: lowest degree polynomials required
    :param: max_degree: highest degree polynomials required
    """

    def __new__(cls, min_degree, max_degree):
        if min_degree == 0:
            # if the restriction is trivial return the original space
            return PolynomialSpace(max_degree, max_degree)
        else:
            return super(RestrictedPolynomialSpace, cls).__new__(cls)

    def __init__(self, min_degree, max_degree):
        self.min_degree = min_degree
        self.max_degree = max_degree

        super(RestrictedPolynomialSpace, self).__init__(-1, max_degree)

    def __repr__(self):
        return "P" + "(min " + str(self.min_degree) + " max " + str(self.max_degree) + ")"

    def to_ON_polynomial_set(self, ref_el):
        if not isinstance(ref_el, reference_element.Cell):
            ref_el = CellComplexToFiat(ref_el)
        sd = ref_el.get_spatial_dimension()

        dimPmin = expansions.polynomial_dimension(ref_el, self.min_degree)
        dimPmax = expansions.polynomial_dimension(ref_el, self.max_degree)

        base_ON = polynomial_set.ONPolynomialSet(ref_el, self.max_degree, (sd,))

        indices = list(chain(*(range(i * dimPmin, i * dimPmax) for i in range(sd))))
        restricted_ON = base_ON.take(indices)
        return restricted_ON



         


class ConstructedPolynomialSpace(PolynomialSpace):
    """
    Sub degree is inherited from the largest of the component spaces,
    super degree is unknown.
    """
    def __init__(self, weights, spaces):

        self.weights = weights
        self.spaces = spaces

        subdegree = max([space.subdegree for space in spaces])

        super(ConstructedPolynomialSpace, self).__init__(subdegree, None)

    def __repr__(self):
        return "+".join([str(w) + "*" + str(x) for (w, x) in zip(self.weights, self.spaces)])

    def to_ON_polynomial_set(self, ref_el):
        if not isinstance(ref_el, reference_element.Cell):
            ref_el = CellComplexToFiat(ref_el)
        space_poly_sets = [s.to_ON_polynomial_set(ref_el) for s in self.spaces]
        sd = ref_el.get_spatial_dimension()

        if all([w == 1 for w in self.weights]):
            return polynomial_set.polynomial_set_union_normalized(space_poly_sets)
        # otherwise have to work on this through tabulation

        k = max([s.superdegree for s in self.spaces])
        Q = create_quadrature(ref_el, 2 * (k + 1))
        Qpts, Qwts = Q.get_points(), Q.get_weights()
        for (space, w) in zip(space_poly_sets, self.weights):
            space_at_Qpts = space.tabulate(Qpts)[(0,) * sd]
            tabulated_expr = tabulate_sympy(w, Qpts).T

            scaled_at_Qpts = space_at_Qpts[:, None, :] * tabulated_expr[None, :, :]
            print(scaled_at_Qpts)

        raise NotImplementedError("constructed not implemented yet")
        return super().to_ON_polynomial_set(ref_el)

    def __mul__(self, x):
        return ConstructedPolynomialSpace([x*w for w in self.weights],
                                          self.spaces)

    def __add__(self, x):
        return ConstructedPolynomialSpace(self.weights.extend([1]),
                                          self.spaces.extend(x))


class VectorPolynomialSpace(PolynomialSpace):

    def __init__(self, *spaces):
        self.component_spaces = []
        for space in spaces:
            assert isinstance(space, PolynomialSpace)
            self.component_spaces.append(space)
            # if isinstance(space, VectorPolynomialSpace):
            #     self.component_spaces.extend(space.component_spaces)

    def dim(self):
        return len(self.component_spaces)

    def complete(self):
        return all([c.complete for c in self.component_spaces])


P0 = PolynomialSpace(0, 0)
P1 = PolynomialSpace(1, 1)
P2 = PolynomialSpace(2, 2)
P3 = PolynomialSpace(3, 3)
P4 = PolynomialSpace(4, 4)

Q1 = PolynomialSpace(1, 2)
Q2 = PolynomialSpace(2, 3)
Q3 = PolynomialSpace(3, 4)
Q4 = PolynomialSpace(4, 5)
