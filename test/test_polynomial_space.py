from redefining_fe import *
import sympy as sp
from redefining_fe.spaces.polynomial_spaces import PolynomialSpace, RestrictedPolynomialSpace, ConstructedPolynomialSpace
from FIAT.polynomial_set import ONPolynomialSet
from FIAT import expansions, polynomial_set, reference_element
from redefining_fe.cells import CellComplexToFiat
from itertools import chain
from FIAT.quadrature_schemes import create_quadrature
import numpy as np


def test_instantiation():
    cell = n_sided_polygon(3)
    x = sp.Symbol("x")
    composite = P0 + x*P0
    on_set = composite.to_ON_polynomial_set(cell)

    assert isinstance(on_set, ONPolynomialSet)


def test_restriction():
    cell = n_sided_polygon(3)
    restricted = P3.restrict(2, 3)

    # doesn't contain constants
    assert restricted.subdegree == -1
    assert restricted.superdegree == 3

    x = sp.Symbol("x")
    composite = P2 + x*restricted
    assert isinstance(composite, ConstructedPolynomialSpace)
    res_on_set = restricted.to_ON_polynomial_set(cell)
    P3_on_set = P3.to_ON_polynomial_set(cell)
    assert res_on_set.get_num_members() < P3_on_set.get_num_members()
    print(composite.spaces)

    not_restricted = P3.restrict(0, 3)
    assert isinstance(not_restricted, PolynomialSpace)
    assert not isinstance(not_restricted, RestrictedPolynomialSpace)


def test_compare_pk_pkp1():
    cell = n_sided_polygon(3)
    ref_el = CellComplexToFiat(cell)
    degree = 2
    sd = ref_el.get_spatial_dimension()

    k = degree - 1
    vec_Pkp1 = polynomial_set.ONPolynomialSet(ref_el, k + 1, (sd,))
    vec_Pk = polynomial_set.ONPolynomialSet(ref_el, k, (sd,))

    dimPkp1 = expansions.polynomial_dimension(ref_el, k + 1)
    dimPk = expansions.polynomial_dimension(ref_el, k)
    dimPkm1 = expansions.polynomial_dimension(ref_el, k - 1)

    vec_Pk_indices = list(chain(*(range(i * dimPkp1, i * dimPkp1 + dimPk)
                                  for i in range(sd))))
    vec_Pk_from_Pkp1 = vec_Pkp1.take(vec_Pk_indices)

    Q = create_quadrature(ref_el, 3)

    def func(x):
        return x[0]
    
    Pkp1 = polynomial_set.ONPolynomialSet(ref_el, k + 1)
    PkH = Pkp1.take(list(range(dimPkm1, dimPk)))

    Q = create_quadrature(ref_el, 2 * (k + 1))
    Qpts, Qwts = Q.get_points(), Q.get_weights()

    # have to work on this through "tabulate" interface
    # first, tabulate PkH at quadrature points
    PkH_at_Qpts = PkH.tabulate(Qpts)[(0,) * sd]
    Pkp1_at_Qpts = Pkp1.tabulate(Qpts)[(0,) * sd]

    x = Qpts.T
    print(PkH_at_Qpts.shape)
    print(x.shape)
    print(PkH_at_Qpts[:, None, :])
    PkHx_at_Qpts = PkH_at_Qpts[:, None, :] * x[None, :, :]
    assert np.allclose(project(func, vec_Pk, Q), project(func, vec_Pk_from_Pkp1, Q))


def project(f, U, Q):
    """Computes the expansion coefficients of f in terms of the members of
    a polynomial set U.  Numerical integration is performed by
    quadrature rule Q.
    """
    pts = Q.get_points()
    wts = Q.get_weights()
    f_at_qps = [f(x) for x in pts]
    U_at_qps = U.tabulate_new(pts)
    coeffs = np.array([sum(wts * f_at_qps * phi) for phi in U_at_qps])
    return coeffs
