from redefining_fe import *
import sympy as sp
from redefining_fe.spaces.polynomial_spaces import PolynomialSpace, RestrictedPolynomialSpace, ConstructedPolynomialSpace
from FIAT.polynomial_set import ONPolynomialSet
from FIAT import expansions, polynomial_set
from redefining_fe.cells import CellComplexToFiat
from itertools import chain
from FIAT.quadrature_schemes import create_quadrature
from redefining_fe.utils import tabulate_sympy, max_deg_sp_mat
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


def test_coeffs():
    cell = n_sided_polygon(3)
    ref_el = CellComplexToFiat(cell)
    sd = ref_el.get_spatial_dimension()
    deg = 2
    k = deg - 1
    vec_Pkp1 = polynomial_set.ONPolynomialSet(ref_el, k + 1, (sd,))
    dimPkp1 = expansions.polynomial_dimension(ref_el, k + 1)
    dimPk = expansions.polynomial_dimension(ref_el, k)

    vec_Pk_indices = list(chain(*(range(i * dimPkp1, i * dimPkp1 + dimPk)
                                  for i in range(sd))))
    vec_Pk_from_Pkp1 = vec_Pkp1.take(vec_Pk_indices)

    vec_Pk = PolynomialSpace(deg - 1, deg - 1, vec=True).to_ON_polynomial_set(ref_el)

    print("A", vec_Pk_from_Pkp1.coeffs)
    print("B", vec_Pk.coeffs)
    from FIAT.polynomial_set import construct_new_coeffs

    print(construct_new_coeffs(ref_el, vec_Pk_from_Pkp1, vec_Pk))


@pytest.mark.parametrize("deg", [1, 2, 3, 4])
def test_scaled_construction(deg):
    cell = n_sided_polygon(3)
    ref_el = CellComplexToFiat(cell)
    sd = ref_el.get_spatial_dimension()
    x = sp.Symbol("x")
    y = sp.Symbol("y")
    M = sp.Matrix([[x, y]])

    vec_Pd = PolynomialSpace(deg - 1, deg - 1, vec=True)
    Pd = PolynomialSpace(deg - 1, deg - 1)
    composite = vec_Pd + (Pd.restrict(deg - 2, deg - 1))*M

    for space in composite.spaces:
        print(space)

    assert isinstance(composite, ConstructedPolynomialSpace)
    on_set = composite.to_ON_polynomial_set(cell)

    from FIAT.raviart_thomas import RTSpace
    rt_space = RTSpace(ref_el, deg)

    Q = create_quadrature(ref_el, 2*(deg+1))
    Qpts, _ = Q.get_points(), Q.get_weights()
    fiat_vals = rt_space.tabulate(Qpts)[(0,) * sd]
    my_vals = on_set.tabulate(Qpts)[(0,) * sd]
    assert np.allclose(fiat_vals, my_vals)
    fiat_vals = svd_reshape(fiat_vals)
    my_vals = svd_reshape(my_vals)

    (x, res, _, _) = np.linalg.lstsq(fiat_vals.T, my_vals.T)
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
    composite = vec_Pk + (Pk.restrict(deg - 2, deg - 1))*M

    assert isinstance(composite, ConstructedPolynomialSpace)

    from FIAT.nedelec import NedelecSpace2D
    fiat_nd_space = NedelecSpace2D(ref_el, deg)
    my_nd_space = composite.to_ON_polynomial_set(cell)

    fiat_Pk, fiat_Pkp1, fiat_crossX, fiat_xP = fiat_nedelec(deg, ref_el)
    my_Pkp1, my_crossX, my_xP = my_nedelec(deg, ref_el, M)
    my_Pk = vec_Pk.to_ON_polynomial_set(ref_el)

    # Creation of space scaled by x matches
    assert (np.allclose(fiat_Pkp1, my_Pkp1))
    assert (np.allclose(fiat_crossX, my_crossX))
    assert (np.allclose(fiat_xP.coeffs, my_xP.coeffs))

    # Creation of Vector Pk space matches
    (_, _, dim) = my_Pk.coeffs.shape
    for i in range(dim):
        assert (np.allclose(fiat_Pk.coeffs[:, :, i], my_Pk.coeffs[:, :, i]))

    # Combination does not - therefore issue is in polynomial_set_union_normalized
    print("Combined Coeffs:", np.allclose(fiat_nd_space.coeffs, my_nd_space.coeffs))

    # The coefficient union is not the problem
    from FIAT.polynomial_set import construct_new_coeffs
    fiat_added_coeffs = construct_new_coeffs(ref_el, fiat_Pk, fiat_xP)
    my_added_coeffs = construct_new_coeffs(ref_el, my_Pk, my_xP)
    assert np.allclose(fiat_added_coeffs, my_added_coeffs)

    # Reshape for svd
    fiat_reshaped = svd_reshape(fiat_added_coeffs)
    my_reshaped = svd_reshape(my_added_coeffs)
    assert np.allclose(fiat_reshaped, my_reshaped)

    # check conditioning
    print(np.linalg.cond(fiat_reshaped))
    print(np.linalg.cond(my_reshaped))
    print(np.allclose(np.linalg.cond(fiat_reshaped), np.linalg.cond(my_reshaped)))

    # Do svd
    fiat_u, fiat_sigma, fiat_svd = np.linalg.svd(fiat_reshaped, 1)
    my_u, my_sigma, my_svd = np.linalg.svd(my_reshaped, 1)
    # fiat_q, _ = np.linalg.qr(fiat_reshaped)
    # my_q, _ = np.linalg.qr(my_reshaped)
    # print(fiat_q - my_q)
    # print(np.allclose(fiat_q, my_q))
    _, n = fiat_u.shape
    m, _ = fiat_svd.shape
    slen = len(fiat_sigma)
    num_sv = len([s for s in my_sigma if abs(s) > 1.e-10])

    slen = len(fiat_sigma)
    fiat_num_sv = len([s for s in fiat_sigma if abs(s) > 1.e-10])
    Q = create_quadrature(ref_el, 2*deg + 2)
    Qpts, _ = Q.get_points(), Q.get_weights()
    print("U:", np.allclose(fiat_u, my_u))
    print("V:", np.allclose(fiat_svd, my_svd, rtol=0, atol=1e-7))
    print('sigmas', np.allclose(fiat_sigma, my_sigma))
    print(fiat_svd - my_svd)
    # print("my v0", my_svd[2])
    # print("fiat v1", fiat_svd[2])
    # print(np.allclose(abs(my_svd[2]), abs(fiat_svd[2])))
    # print("fiat v0", fiat_svd[0])
    # print("my v1", my_svd[1])

    print(np.allclose(fiat_reshaped, np.dot(fiat_u * fiat_sigma, fiat_svd[:fiat_num_sv, :])))
    print(np.allclose(my_reshaped, np.dot(my_u * my_sigma, my_svd[:num_sv, :])))

    # print("test", np.allclose(fiat_nd_space.expansion_set._tabulate(fiat_nd_space.embedded_degree, Qpts, order=0)[(0,) * sd], fiat_u*fiat_sigma))
    my_recon = np.dot(my_u * my_sigma, my_svd[:slen, :])
    fiat_recon = np.dot(fiat_u * fiat_sigma, fiat_svd[:slen, :])
    print(np.allclose(my_recon, fiat_recon))

    # Round then do svd
    _, fiat_sigma, fiat_svd = np.linalg.svd(fiat_reshaped.round(14), 1)
    _, my_sigma, my_svd = np.linalg.svd(my_reshaped.round(14), 1)
    print("SVD (after rounding):", np.allclose(fiat_svd, my_svd))
    print('sigmas', np.allclose(fiat_sigma, my_sigma))

    fiat_vals = fiat_nd_space.tabulate(Qpts)[(0,) * sd]
    my_vals = my_nd_space.tabulate(Qpts)[(0,) * sd]
    # print((fiat_vals - my_vals))
    for r in fiat_nd_space.coeffs:
        print(r in my_nd_space.coeffs)
    print("tabulated", np.allclose(fiat_vals, my_vals))
    print(fiat_vals.shape)
    fiat_vals = svd_reshape(fiat_vals)
    my_vals = svd_reshape(my_vals)
    print(fiat_vals.shape)
    print(my_vals.shape)
    # confirm the two tabulations have the same span
    (x, res, rank, _) = np.linalg.lstsq(fiat_vals.T, my_vals.T)
    (x2, res2, rank, _) = np.linalg.lstsq(my_vals.T, fiat_vals.T)
    # print(x.shape)
    x1 = np.linalg.inv(x)
    print(np.linalg.norm(x1))
    assert np.allclose(np.linalg.norm(my_vals.T - fiat_vals.T @ x), 0)
    assert np.allclose(np.linalg.norm(fiat_vals.T - my_vals.T @ x1), 0)
    assert np.allclose(res, 0)
    assert np.allclose(res2, 0)
    print(res2)
    print(res)
    print(rank)
    print(rank)


def svd_reshape(new_coeffs):
    func_shape = new_coeffs.shape[1:]
    new_shape0 = new_coeffs.shape[0]
    new_shape1 = np.prod(func_shape)
    newshape = (new_shape0, new_shape1)
    nc = np.reshape(new_coeffs, newshape)
    return nc


def fiat_nedelec(deg, ref_el):
    sd = ref_el.get_spatial_dimension()
    k = deg - 1
    Q = create_quadrature(ref_el, 2 * (k + 1))
    Qpts, Qwts = Q.get_points(), Q.get_weights()

    vec_Pkp1 = polynomial_set.ONPolynomialSet(ref_el, k + 1, (sd,), scale="orthonormal")
    dimPkp1 = expansions.polynomial_dimension(ref_el, k + 1)
    dimPk = expansions.polynomial_dimension(ref_el, k)
    dimPkm1 = expansions.polynomial_dimension(ref_el, k - 1)

    vec_Pk_indices = list(chain(*(range(i * dimPkp1, i * dimPkp1 + dimPk)
                                  for i in range(sd))))
    vec_Pk_from_Pkp1 = vec_Pkp1.take(vec_Pk_indices)

    Pkp1 = polynomial_set.ONPolynomialSet(ref_el, k + 1, scale="orthonormal")
    PkH = Pkp1.take(list(range(dimPkm1, dimPk)))

    PkH_at_Qpts = PkH.tabulate(Qpts)[(0,) * sd]
    Pkp1_at_Qpts = Pkp1.tabulate(Qpts)[(0,) * sd]

    CrossX = np.dot(np.array([[0.0, 1.0], [-1.0, 0.0]]), Qpts.T)
    PkHCrossX_at_Qpts = PkH_at_Qpts[:, None, :] * CrossX[None, :, :]
    PkHCrossX_coeffs = np.dot(np.multiply(PkHCrossX_at_Qpts, Qwts), Pkp1_at_Qpts.T)
    PkHcrossX = polynomial_set.PolynomialSet(ref_el,
                                             k + 1,
                                             k + 1,
                                             vec_Pkp1.get_expansion_set(),
                                             PkHCrossX_coeffs)
    return vec_Pk_from_Pkp1, Pkp1_at_Qpts, PkHCrossX_coeffs, PkHcrossX


def my_nedelec(deg, ref_el, M):
    sd = ref_el.get_spatial_dimension()
    Pd = PolynomialSpace(deg - 1, deg - 1)
    w_deg = max_deg_sp_mat(M)
    space = (Pd.restrict(deg - 2, deg - 1)).to_ON_polynomial_set(ref_el)
    Pkpw = polynomial_set.ONPolynomialSet(ref_el, space.degree + w_deg, scale="orthonormal")
    vec_Pkpw = polynomial_set.ONPolynomialSet(ref_el, space.degree + w_deg, (sd,), scale="orthonormal")

    k = deg - 1
    Q = create_quadrature(ref_el, 2 * (k + 1))
    Qpts, Qwts = Q.get_points(), Q.get_weights()
    space_at_Qpts = space.tabulate(Qpts)[(0,) * sd]
    Pkpw_at_Qpts = Pkpw.tabulate(Qpts)[(0,) * sd]

    tabulated_expr = tabulate_sympy(M, Qpts).T
    scaled_at_Qpts = space_at_Qpts[:, None, :] * tabulated_expr[None, :, :]
    PkHw_coeffs = np.dot(np.multiply(scaled_at_Qpts, Qwts), Pkpw_at_Qpts.T)
    on = polynomial_set.PolynomialSet(ref_el, space.degree + w_deg,
                                      space.degree + w_deg,
                                      vec_Pkpw.get_expansion_set(),
                                      PkHw_coeffs)
    return Pkpw_at_Qpts, PkHw_coeffs, on


def test_embedding():
    cell = n_sided_polygon(3)
    ref_el = CellComplexToFiat(cell)
    sd = ref_el.get_spatial_dimension()

    A = P2.to_ON_polynomial_set(cell)
    B = P3.to_ON_polynomial_set(cell)
    print(sd)

    print(A.degree)
    print(B.degree)

    new_coeffs = polynomial_set.construct_new_coeffs(ref_el, A, B)
    print(new_coeffs)


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
