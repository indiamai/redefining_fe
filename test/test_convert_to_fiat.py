import pytest
import numpy as np
# import sympy as sp
from redefining_fe import *
from firedrake import *
from FIAT.quadrature_schemes import create_quadrature
from ufl.cell import simplex
from test_2d_examples_docs import construct_nd, construct_rt, construct_cg3
from test_polynomial_space import flatten

vert = Point(0)
edge = Point(1, [Point(0), Point(0)], vertex_num=2)
tri = n_sided_polygon(3)


def create_dg1(cell):
    xs = [DOF(DeltaPairing(), PointKernel(cell.vertices(return_coords=True)[0]))]
    Pk = PolynomialSpace(1)
    dg = ElementTriple(cell, (Pk, CellL2, C0), DOFGenerator(xs, get_cyc_group(len(cell.vertices())), S1))
    return dg


def create_dg2(cell):
    xs = [DOF(DeltaPairing(), PointKernel(cell.vertices(return_coords=True)[0]))]
    center = [DOF(DeltaPairing(), PointKernel((0,)))]
    Pk = PolynomialSpace(2)
    dg = ElementTriple(cell, (Pk, CellL2, C0), [DOFGenerator(xs, get_cyc_group(len(cell.vertices())), S1),
                                                DOFGenerator(center, S1, S1)])
    return dg


def create_dg1_uneven(cell):
    xs = [DOF(DeltaPairing(), PointKernel(-0.75,))]
    center = [DOF(DeltaPairing(), PointKernel((0.25,)))]
    Pk = PolynomialSpace(1)
    dg = ElementTriple(cell, (Pk, CellL2, C0), [DOFGenerator(xs, S1, S2),
                                                DOFGenerator(center, S1, S2)])
    return dg


def create_cg1(cell):
    deg = 1
    vert_dg = create_dg1(cell.vertices(get_class=True)[0])
    xs = [immerse(cell, vert_dg, TrH1)]

    Pk = PolynomialSpace(deg)
    cg = ElementTriple(cell, (Pk, CellL2, C0), DOFGenerator(xs, get_cyc_group(len(cell.vertices())), S1))
    return cg


def create_cg1_flipped(cell):
    deg = 1
    vert_dg = create_dg1(cell.vertices(get_class=True)[0])
    xs = [immerse(cell, vert_dg, TrH1, node=1)]

    Pk = PolynomialSpace(deg)
    cg = ElementTriple(cell, (Pk, CellL2, C0), DOFGenerator(xs, get_cyc_group(len(cell.vertices())), S1))

    for dof in cg.generate():
        print(dof)
    return cg


def create_cg2(cell):
    deg = 2
    if cell.dim() > 1:
        raise NotImplementedError("This method is for cg2 on edges, please use create_cg2_tri for triangles")
    vert_dg = create_dg1(cell.vertices(get_class=True)[0])
    xs = [immerse(cell, vert_dg, TrH1)]
    center = [DOF(DeltaPairing(), PointKernel((0,)))]

    Pk = PolynomialSpace(deg)
    cg = ElementTriple(cell, (Pk, CellL2, C0), [DOFGenerator(xs, get_cyc_group(len(cell.vertices())), S1),
                                                DOFGenerator(center, S1, S1)])
    return cg


def create_cg2_tri(cell):
    deg = 2
    Pk = PolynomialSpace(deg)

    vert_dg0 = create_dg1(cell.vertices(get_class=True)[0])
    xs = [immerse(cell, vert_dg0, TrH1)]

    edge_dg0 = ElementTriple(cell.edges(get_class=True)[0], (Pk, CellL2, C0), DOFGenerator([DOF(DeltaPairing(), PointKernel((0,)))], S1, S1))
    edge_xs = [immerse(cell, edge_dg0, TrH1)]

    cg = ElementTriple(cell, (Pk, CellL2, C0), [DOFGenerator(xs, get_cyc_group(len(cell.vertices())), S1),
                                                DOFGenerator(edge_xs, tri_C3, S1)])
    return cg


@pytest.mark.parametrize("cell", [tri])
def test_create_fiat_nd(cell):
    nd = construct_nd(cell)
    ref_el = cell.to_fiat()
    sd = ref_el.get_spatial_dimension()
    deg = 1

    from FIAT.nedelec import Nedelec
    fiat_elem = Nedelec(ref_el, deg)
    my_elem = nd.to_fiat_elem()

    Q = create_quadrature(ref_el, 2*(deg+1))
    Qpts, _ = Q.get_points(), Q.get_weights()

    fiat_vals = fiat_elem.tabulate(0, Qpts)
    my_vals = my_elem.tabulate(0, Qpts)

    fiat_vals = flatten(fiat_vals[(0,) * sd])
    my_vals = flatten(my_vals[(0,) * sd])

    (x, res, _, _) = np.linalg.lstsq(fiat_vals.T, my_vals.T)
    x1 = np.linalg.inv(x)
    assert np.allclose(np.linalg.norm(my_vals.T - fiat_vals.T @ x), 0)
    assert np.allclose(np.linalg.norm(fiat_vals.T - my_vals.T @ x1), 0)
    assert np.allclose(res, 0)


@pytest.mark.parametrize("cell", [tri])
def test_create_fiat_rt(cell):
    rt = construct_rt(cell)
    ref_el = cell.to_fiat()
    sd = ref_el.get_spatial_dimension()
    deg = 1

    from FIAT.raviart_thomas import RaviartThomas
    fiat_elem = RaviartThomas(ref_el, deg)
    my_elem = rt.to_fiat_elem()

    Q = create_quadrature(ref_el, 2*(deg+1))
    Qpts, _ = Q.get_points(), Q.get_weights()

    fiat_vals = fiat_elem.tabulate(0, Qpts)
    my_vals = my_elem.tabulate(0, Qpts)

    fiat_vals = flatten(fiat_vals[(0,) * sd])
    my_vals = flatten(my_vals[(0,) * sd])

    (x, res, _, _) = np.linalg.lstsq(fiat_vals.T, my_vals.T)
    x1 = np.linalg.inv(x)
    assert np.allclose(np.linalg.norm(my_vals.T - fiat_vals.T @ x), 0)
    assert np.allclose(np.linalg.norm(fiat_vals.T - my_vals.T @ x1), 0)
    assert np.allclose(res, 0)


@pytest.mark.parametrize("elem_gen,elem_code,deg,cell", [(create_cg1, "CG", 1, edge),
                                                         (create_dg1, "DG", 1, edge),
                                                         (create_dg2, "DG", 2, edge),
                                                         (create_cg2, "CG", 2, edge)])
def test_create_fiat_lagrange(elem_gen, elem_code, deg, cell):
    elem = elem_gen(cell)
    ref_el = cell.to_fiat()
    sd = ref_el.get_spatial_dimension()

    from FIAT.lagrange import Lagrange
    fiat_elem = Lagrange(ref_el, deg)

    my_elem = elem.to_fiat_elem()

    Q = create_quadrature(ref_el, 2*(deg+1))
    Qpts, _ = Q.get_points(), Q.get_weights()

    fiat_vals = fiat_elem.tabulate(0, Qpts)
    my_vals = my_elem.tabulate(0, Qpts)

    fiat_vals = flatten(fiat_vals[(0,) * sd])
    my_vals = flatten(my_vals[(0,) * sd])

    (x, res, _, _) = np.linalg.lstsq(fiat_vals.T, my_vals.T)
    x1 = np.linalg.inv(x)
    assert np.allclose(np.linalg.norm(my_vals.T - fiat_vals.T @ x), 0)
    assert np.allclose(np.linalg.norm(fiat_vals.T - my_vals.T @ x1), 0)
    assert np.allclose(res, 0)


@pytest.mark.parametrize("elem_gen,elem_code,deg,cell", [(create_cg1, "CG", 1, edge),
                                                         (create_dg1, "DG", 1, edge),
                                                         (create_dg2, "DG", 2, edge),
                                                         (create_cg2, "CG", 2, edge),
                                                         (create_cg2_tri, "CG", 2, tri),
                                                         (create_cg1, "CG", 1, tri),
                                                         (create_dg1, "DG", 1, tri),
                                                         (construct_cg3, "CG", 3, tri),
                                                         (construct_nd, "N1curl", 1, tri)])
def test_entity_perms(elem_gen, elem_code, deg, cell):
    elem = elem_gen(cell)

    print(elem.to_fiat_elem())


@pytest.mark.parametrize("elem_gen,elem_code,deg", [(create_cg1, "CG", 1),
                                                    (create_dg1, "DG", 1),
                                                    (create_dg2, "DG", 2),
                                                    (create_cg2, "CG", 2)
                                                    ])
def test_2d(elem_gen, elem_code, deg):
    cell = edge
    elem = elem_gen(cell)

    mesh = UnitIntervalMesh(5)
    V1 = FunctionSpace(mesh, elem_code, deg)
    u = TrialFunction(V1)
    v = TestFunction(V1)
    f = Function(V1)
    x, = SpatialCoordinate(mesh)
    f.interpolate(cos(x*pi*2))
    a = (inner(grad(u), grad(v)) + inner(u, v)) * dx
    L = inner(f, v) * dx
    u1 = Function(V1)
    solve(a == L, u1)

    V2 = FunctionSpace(mesh, elem.to_ufl_elem())
    u = TrialFunction(V2)
    v = TestFunction(V2)
    f = Function(V2)
    x, = SpatialCoordinate(mesh)
    f.interpolate((1+8*pi*pi)*cos(x*pi*2))
    a = (inner(grad(u), grad(v)) + inner(u, v)) * dx
    L = inner(f, v) * dx
    u2 = Function(V2)
    solve(a == L, u2)

    res = sqrt(assemble(dot(u1 - u1, u1 - u2) * dx))
    assert np.allclose(res, 0)


# pytest.param( marks=pytest.mark.xfail(reason='Orientations bug'))
# @pytest.mark.parametrize("elem_gen,elem_code,deg", [pytest.param(create_cg1, "CG", 1, marks=pytest.mark.xfail(reason='Values incorrect')),
#                                                    pytest.param(create_dg1, "DG", 1, marks=pytest.mark.xfail(reason='Values incorrect')),
#                                                    pytest.param(construct_cg3, "CG", 3, marks=pytest.mark.xfail(reason='Firedrake error'))])
@pytest.mark.parametrize("elem_gen,elem_code,deg,conv_rate", [(create_cg1, "CG", 1, 1.8)])
def test_helmholtz(elem_gen, elem_code, deg, conv_rate):
    cell = n_sided_polygon(3)
    elem = elem_gen(cell)

    diff_fiat = [0 for i in range(3, 6)]
    diff_new = [0 for i in range(3, 6)]
    for i in range(3, 6):
        mesh = UnitSquareMesh(2 ** i, 2 ** i)

        V = FunctionSpace(mesh, elem_code, deg)
        res1 = helmholtz_solve(mesh, V)
        diff_fiat[i - 3] = res1

        V = FunctionSpace(mesh, elem.to_ufl_elem())
        res2 = helmholtz_solve(mesh, V)
        diff_new[i - 3] = res2

    diff_fiat = np.array(diff_fiat)
    print("l2 error norms:", diff_fiat)
    conv = np.log2(diff_fiat[:-1] / diff_fiat[1:])
    print("convergence order:", conv)
    assert (np.array(conv) > conv_rate).all()

    print("l2 error norms:", diff_new)
    diff_new = np.array(diff_new)
    conv = np.log2(diff_new[:-1] / diff_new[1:])
    print("convergence order:", conv)
    assert (np.array(conv) > conv_rate).all()


def helmholtz_solve(mesh, V):
    u = TrialFunction(V)
    v = TestFunction(V)
    f = Function(V)
    x, y = SpatialCoordinate(mesh)
    f.interpolate((1+8*pi*pi)*cos(x*pi*2)*cos(y*pi*2))
    a = (inner(grad(u), grad(v)) + inner(u, v)) * dx
    L = inner(f, v) * dx
    u = Function(V)
    solve(a == L, u)
    f.interpolate(cos(x*pi*2)*cos(y*pi*2))
    return sqrt(assemble(dot(u - f, u - f) * dx))


@pytest.mark.parametrize("cell", [vert, edge, tri])
def test_ufl_cell_conversion(cell):
    existing_cell = simplex(len(cell.vertices()))
    print(type(existing_cell))
    ufl_cell = cell.to_ufl()
    print(isinstance(ufl_cell, ufl.Cell))
    print(ufl_cell.cell_complex)
    print(ufl_cell.cellname())


@pytest.mark.parametrize("cell", [edge])
def test_functional_evaluation(cell):
    cg = create_cg1(cell)
    cg_f = create_cg1_flipped(cell)
    ref_el = cell.to_fiat()
    deg = 1

    from FIAT.lagrange import Lagrange
    fiat_elem = Lagrange(ref_el, deg)
    my_elem = cg.to_fiat_elem()
    my_elem_f = cg_f.to_fiat_elem()

    print([n.pt_dict for n in my_elem.dual.nodes])
    print([n.pt_dict for n in my_elem_f.dual.nodes])
    print([n.pt_dict for n in fiat_elem.dual.nodes])

    print("my poly set")
    print(np.matmul(my_elem.V, my_elem.get_coeffs().T))
    print(np.matmul(my_elem_f.V, my_elem.get_coeffs().T))
    # print(np.matmul(fiat_elem.V.T, my_elem.get_coeffs()))

    print("my poly set")
    print(np.matmul(my_elem.V, my_elem_f.get_coeffs().T))
    print(np.matmul(my_elem_f.V, my_elem_f.get_coeffs().T))


@pytest.mark.parametrize("cell", [edge])
def test_functional_evaluation_uneven(cell):
    cg = create_dg1(cell)
    cg_f = create_dg1_uneven(cell)

    print("EVEN")
    my_elem = cg.to_fiat_elem()
    print("UNEVEN")
    my_elem_f = cg_f.to_fiat_elem()
    print(my_elem_f)
    print(my_elem)


# @pytest.mark.parametrize("cell", [pytest.param(tri, marks=pytest.mark.xfail(reason="Dense matrix dimensions in vector case"))])
@pytest.mark.parametrize("cell", [tri])
def test_functional_evaluation_vector(cell):
    rt = construct_rt(cell)

    from FIAT.raviart_thomas import RaviartThomas
    ref_el = cell.to_fiat()
    deg = 1
    fiat_elem = RaviartThomas(ref_el, deg)
    my_elem = rt.to_fiat_elem()
    print(my_elem)
    print(fiat_elem)
    # deg = 1

    # x = sp.Symbol("x")
    # y = sp.Symbol("y")

    # M = sp.Matrix([[x, y]])
    # vec_Pd = PolynomialSpace(deg - 1, set_shape=True)
    # Pd = PolynomialSpace(deg - 1)
    # rt_space = vec_Pd + (Pd.restrict(deg - 2, deg - 1))*M

    # tri = n_sided_polygon(3)
    # edge = tri.edges(get_class=True)[0]

    # xs = [DOF(L2InnerProd(), PolynomialKernel(1))]
    # dofs = DOFGenerator(xs, S1, S2)

    # int_rt = ElementTriple(edge, (P1, CellHDiv, C0), dofs)

    # int_rt.to_fiat_elem()
