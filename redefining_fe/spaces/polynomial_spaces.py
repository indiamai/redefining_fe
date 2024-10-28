from FIAT.polynomial_set import ONPolynomialSet
from FIAT.quadrature_schemes import create_quadrature
from FIAT import expansions, polynomial_set, reference_element
from itertools import chain
from redefining_fe.utils import tabulate_sympy, max_deg_sp_mat
import sympy as sp
import numpy as np


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

    def __init__(self, subdegree, superdegree, vec=False):
        self.subdegree = subdegree
        self.superdegree = superdegree
        self.vec = vec

    def complete(self):
        return self.subdegree == self.superdegree

    def degree(self):
        return self.subdegree

    def to_ON_polynomial_set(self, ref_el, k=None):
        # how does super/sub degrees work here
        if not isinstance(ref_el, reference_element.Cell):
            ref_el = ref_el.to_fiat()
        sd = ref_el.get_spatial_dimension()
        if self.vec:
            shape = (sd,)
        else:
            shape = tuple()

        return ONPolynomialSet(ref_el, self.subdegree, shape, scale="orthonormal")

    def __repr__(self):
        res = ""
        if self.complete():
            res += "P" + str(self.subdegree)
        else:
            res += "Psub" + str(self.subdegree) + "sup" + str(self.superdegree)
        if self.vec:
            res += "^d"
        return res

    def __mul__(self, x):
        """
        When multiplying a Polynomial Space by a sympy object, you need to multiply with
        the sympy object on the right. This is due to Sympy's implementation of __mul__ not
        passing to this handler as it should.
        """
        if isinstance(x, sp.Symbol):
            return ConstructedPolynomialSpace([x], [self])
        elif isinstance(x, sp.Matrix):
            return ConstructedPolynomialSpace([x], [self])
        else:
            raise TypeError(f'Cannot multiply a PolySpace with {type(x)}')

    __rmul__ = __mul__

    def __add__(self, x):
        return ConstructedPolynomialSpace([1, 1], [self, x])

    def restrict(self, min_degree, max_degree):
        return RestrictedPolynomialSpace(min_degree, max_degree, self.vec)

    def _to_dict(self):
        return {"vec": self.vec, "sub": self.subdegree, "super": self.superdegree}

    def dict_id(self):
        return "PolynomialSpace"

    def _from_dict(obj_dict):
        return PolynomialSpace(obj_dict["sub"], obj_dict["super"], obj_dict["vec"])


class RestrictedPolynomialSpace(PolynomialSpace):
    """
    Represents a polynomial space of all polynomials between two degrees.

    :param: min_degree: lowest degree polynomials required (-1 to include constants)
    :param: max_degree: highest degree polynomials required
    """

    def __new__(cls, min_degree, max_degree, vec=False):
        if min_degree == -1:
            # if the restriction is trivial return the original space
            return PolynomialSpace(max_degree, max_degree, vec)
        else:
            return super(RestrictedPolynomialSpace, cls).__new__(cls)

    def __init__(self, min_degree, max_degree, vec=False):
        self.min_degree = min_degree
        self.max_degree = max_degree

        super(RestrictedPolynomialSpace, self).__init__(-1, max_degree, vec)

    def __repr__(self):
        res = "P" + "(min " + str(self.min_degree) + " max " + str(self.max_degree) + ")"
        if self.vec:
            res += "^d"
        return res

    def to_ON_polynomial_set(self, ref_el):
        if not isinstance(ref_el, reference_element.Cell):
            ref_el = ref_el.to_fiat()
        sd = ref_el.get_spatial_dimension()

        dimPmin = expansions.polynomial_dimension(ref_el, self.min_degree)
        dimPmax = expansions.polynomial_dimension(ref_el, self.max_degree)

        if self.vec:
            base_ON = polynomial_set.ONPolynomialSet(ref_el, self.max_degree, (sd,), scale="orthonormal")
        else:
            base_ON = polynomial_set.ONPolynomialSet(ref_el, self.max_degree, scale="orthonormal")

        indices = list(chain(*(range(i * dimPmin, i * dimPmax) for i in range(sd))))
        restricted_ON = base_ON.take(indices)
        return restricted_ON

    def _to_dict(self):
        super_dict = super(RestrictedPolynomialSpace, self)._to_dict()
        super_dict["min_degree"] = self.min_degree
        super_dict["max_degree"] = self.max_degree
        return super_dict

    def dict_id(self):
        return "RestrictedPolynomialSpace"

    def _from_dict(obj_dict):
        return RestrictedPolynomialSpace(obj_dict["min_degree"], obj_dict["max_degree"], obj_dict["vec"])


class ConstructedPolynomialSpace(PolynomialSpace):
    """
    Sub degree is inherited from the largest of the component spaces,
    super degree is unknown.

    weights can either be 1 or a polynomial in x, where x in R^d
    """
    def __init__(self, weights, spaces):

        self.weights = weights
        self.spaces = spaces

        subdegree = max([space.subdegree for space in spaces])
        vec = any([s.vec for s in spaces])

        super(ConstructedPolynomialSpace, self).__init__(subdegree, -1, vec=vec)

    def __repr__(self):
        return "+".join([str(w) + "*" + str(x) for (w, x) in zip(self.weights, self.spaces)])

    def to_ON_polynomial_set(self, ref_el):
        if not isinstance(ref_el, reference_element.Cell):
            ref_el = ref_el.to_fiat()
        k = max([s.superdegree for s in self.spaces])
        space_poly_sets = [s.to_ON_polynomial_set(ref_el) for s in self.spaces]
        sd = ref_el.get_spatial_dimension()

        if all([w == 1 for w in self.weights]):
            weighted_sets = space_poly_sets

        # otherwise have to work on this through tabulation

        Q = create_quadrature(ref_el, 2 * (k + 1))
        Qpts, Qwts = Q.get_points(), Q.get_weights()
        weighted_sets = []

        for (space, w) in zip(space_poly_sets, self.weights):
            if not (isinstance(w, sp.Expr) or isinstance(w, sp.Matrix)):
                weighted_sets.append(space)
            else:
                w_deg = max_deg_sp_mat(w)
                Pkpw = polynomial_set.ONPolynomialSet(ref_el, space.degree + w_deg, scale="orthonormal")
                vec_Pkpw = polynomial_set.ONPolynomialSet(ref_el, space.degree + w_deg, (sd,), scale="orthonormal")

                space_at_Qpts = space.tabulate(Qpts)[(0,) * sd]
                Pkpw_at_Qpts = Pkpw.tabulate(Qpts)[(0,) * sd]

                tabulated_expr = tabulate_sympy(w, Qpts).T
                scaled_at_Qpts = space_at_Qpts[:, None, :] * tabulated_expr[None, :, :]
                PkHw_coeffs = np.dot(np.multiply(scaled_at_Qpts, Qwts), Pkpw_at_Qpts.T)
                weighted_sets.append(polynomial_set.PolynomialSet(ref_el,
                                                                  space.degree + w_deg,
                                                                  space.degree + w_deg,
                                                                  vec_Pkpw.get_expansion_set(),
                                                                  PkHw_coeffs))
        combined_sets = weighted_sets[0]
        for i in range(1, len(weighted_sets)):
            combined_sets = polynomial_set.polynomial_set_union_normalized(combined_sets, weighted_sets[i])
        return combined_sets

    def __mul__(self, x):
        return ConstructedPolynomialSpace([x*w for w in self.weights],
                                          self.spaces)
    __rmul__ = __mul__

    def __add__(self, x):
        return ConstructedPolynomialSpace(self.weights.extend([1]),
                                          self.spaces.extend(x))

    def _to_dict(self):
        super_dict = super(ConstructedPolynomialSpace, self)._to_dict()
        super_dict["spaces"] = self.spaces
        super_dict["weights"] = self.weights
        return super_dict

    def dict_id(self):
        return "ConstructedPolynomialSpace"

    def _from_dict(obj_dict):
        return ConstructedPolynomialSpace(obj_dict["weights"], obj_dict["spaces"])


P0 = PolynomialSpace(0, 0)
P1 = PolynomialSpace(1, 1)
P2 = PolynomialSpace(2, 2)
P3 = PolynomialSpace(3, 3)
P4 = PolynomialSpace(4, 4)

Q1 = PolynomialSpace(1, 2)
Q2 = PolynomialSpace(2, 3)
Q3 = PolynomialSpace(3, 4)
Q4 = PolynomialSpace(4, 5)
