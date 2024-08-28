from redefining_fe import *
import sympy as sp
from redefining_fe.spaces.polynomial_spaces import PolynomialSpace, RestrictedPolynomialSpace, ConstructedPolynomialSpace
from FIAT.polynomial_set import ONPolynomialSet
from FIAT import polynomial_set
from redefining_fe.cells import CellComplexToFiat
from FIAT.quadrature_schemes import create_quadrature
import numpy as np
import pytest


def test_instantiation():
    cell = n_sided_polygon(3)

    on_set = P2.to_ON_polynomial_set(cell)
    assert isinstance(on_set, ONPolynomialSet)


def test_unscaled_construction():
    cell = n_sided_polygon(3)
    composite = P0 + P1
    on_set = composite.to_ON_polynomial_set(cell)
    assert isinstance(on_set, polynomial_set.PolynomialSet)

    vec_P0 = PolynomialSpace(0, 0, vec=True)
    vec_P1 = PolynomialSpace(1, 1, vec=True)

    composite = vec_P0 + vec_P1
    on_set = composite.to_ON_polynomial_set(cell)
    assert isinstance(on_set, polynomial_set.PolynomialSet)


def test_restriction():
    cell = n_sided_polygon(3)
    restricted = P3.restrict(2, 3)

    # doesn't contain constants
    assert restricted.subdegree == -1
    assert restricted.superdegree == 3

    res_on_set = restricted.to_ON_polynomial_set(cell)
    P3_on_set = P3.to_ON_polynomial_set(cell)
    assert res_on_set.get_num_members() < P3_on_set.get_num_members()

    not_restricted = P3.restrict(-1, 3)
    assert isinstance(not_restricted, PolynomialSpace)
    assert not isinstance(not_restricted, RestrictedPolynomialSpace)


@pytest.mark.parametrize("deg", [1, 2, 3, 4])
def test_rt_construction(deg):
    cell = n_sided_polygon(3)
    ref_el = CellComplexToFiat(cell)
    sd = ref_el.get_spatial_dimension()
    x = sp.Symbol("x")
    y = sp.Symbol("y")
    M = sp.Matrix([[x, y]])

    vec_Pd = PolynomialSpace(deg - 1, deg - 1, vec=True)
    Pd = PolynomialSpace(deg - 1, deg - 1)
    composite = vec_Pd + (Pd.restrict(deg - 2, deg - 1))*M

    assert isinstance(composite, ConstructedPolynomialSpace)
    on_set = composite.to_ON_polynomial_set(cell)

    from FIAT.raviart_thomas import RTSpace
    rt_space = RTSpace(ref_el, deg)

    Q = create_quadrature(ref_el, 2*(deg+1))
    Qpts, _ = Q.get_points(), Q.get_weights()
    fiat_vals = rt_space.tabulate(Qpts)[(0,) * sd]
    my_vals = on_set.tabulate(Qpts)[(0,) * sd]
    assert np.allclose(fiat_vals, my_vals)
    fiat_vals = flatten(fiat_vals)
    my_vals = flatten(my_vals)

    (x, res, _, _) = np.linalg.lstsq(fiat_vals.T, my_vals.T)
    x1 = np.linalg.inv(x)
    assert np.allclose(np.linalg.norm(my_vals.T - fiat_vals.T @ x), 0)
    assert np.allclose(np.linalg.norm(fiat_vals.T - my_vals.T @ x1), 0)
    assert np.allclose(res, 0)


@pytest.mark.parametrize("deg", [1, 2, 3, 4])
def test_nedelec_construction(deg):
    cell = n_sided_polygon(3)
    ref_el = CellComplexToFiat(cell)
    sd = ref_el.get_spatial_dimension()

    x = sp.Symbol("x")
    y = sp.Symbol("y")
    M = sp.Matrix([[y, -x]])

    vec_Pk = PolynomialSpace(deg - 1, deg - 1, vec=True)
    Pk = PolynomialSpace(deg - 1, deg - 1)
    nd = vec_Pk + (Pk.restrict(deg - 2, deg - 1))*M
    assert isinstance(nd, ConstructedPolynomialSpace)

    from FIAT.nedelec import NedelecSpace2D
    fiat_nd_space = NedelecSpace2D(ref_el, deg)
    my_nd_space = nd.to_ON_polynomial_set(cell)

    Q = create_quadrature(ref_el, 2*deg + 2)
    Qpts, _ = Q.get_points(), Q.get_weights()
    fiat_vals = fiat_nd_space.tabulate(Qpts)[(0,) * sd]
    my_vals = my_nd_space.tabulate(Qpts)[(0,) * sd]
    fiat_vals = flatten(fiat_vals)
    my_vals = flatten(my_vals)

    # confirm the two tabulations have the same span
    (x, res, _, _) = np.linalg.lstsq(fiat_vals.T, my_vals.T)
    x1 = np.linalg.inv(x)
    assert np.allclose(np.linalg.norm(my_vals.T - fiat_vals.T @ x), 0)
    assert np.allclose(np.linalg.norm(fiat_vals.T - my_vals.T @ x1), 0)
    assert np.allclose(res, 0)


def flatten(coeffs):
    func_shape = coeffs.shape[1:]
    new_shape0 = coeffs.shape[0]
    new_shape1 = np.prod(func_shape)
    newshape = (new_shape0, new_shape1)
    nc = np.reshape(coeffs, newshape)
    return nc
