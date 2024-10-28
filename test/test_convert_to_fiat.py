import pytest
import numpy as np
from redefining_fe import *
from FIAT.quadrature_schemes import create_quadrature
from firedrake import *
from ufl.cell import simplex
# from test_2d_examples_docs import construct_cg3

vert = Point(0)
edge = Point(1, [Point(0), Point(0)], vertex_num=2)
tri = n_sided_polygon(3)


def create_dg1(cell):
    xs = [DOF(DeltaPairing(), PointKernel(cell.vertices(return_coords=True)[0]))]
    Pk = PolynomialSpace(1, 1)
    dg = ElementTriple(cell, (Pk, CellL2, C0), DOFGenerator(xs, get_cyc_group(len(cell.vertices())), S1))
    return dg


def create_dg2(cell):
    xs = [DOF(DeltaPairing(), PointKernel(cell.vertices(return_coords=True)[0]))]
    center = [DOF(DeltaPairing(), PointKernel((0,)))]
    Pk = PolynomialSpace(2, 2)
    dg = ElementTriple(cell, (Pk, CellL2, C0), [DOFGenerator(xs, get_cyc_group(len(cell.vertices())), S1),
                                                DOFGenerator(center, S1, S1)])
    return dg


def create_cg1(cell):
    deg = 1
    vert_dg = create_dg1(cell.vertices(get_class=True)[0])
    xs = [immerse(cell, vert_dg, TrH1)]

    Pk = PolynomialSpace(deg, deg)
    cg = ElementTriple(cell, (Pk, CellL2, C0), DOFGenerator(xs, get_cyc_group(len(cell.vertices())), S1))
    return cg


def create_cg2(cell):
    deg = 2
    vert_dg = create_dg1(cell.vertices(get_class=True)[0])
    xs = [immerse(cell, vert_dg, TrH1)]
    center = [DOF(DeltaPairing(), PointKernel((0,)))]

    Pk = PolynomialSpace(deg, deg)
    cg = ElementTriple(cell, (Pk, CellL2, C0), [DOFGenerator(xs, get_cyc_group(len(cell.vertices())), S1),
                                                DOFGenerator(center, S1, S1)])
    return cg


@pytest.mark.parametrize("cell", [tri, edge])
def test_create_fiat_cg1(cell):
    deg = 1
    cg = create_cg1(cell)
    ref_el = cell.to_fiat()
    sd = ref_el.get_spatial_dimension()

    from FIAT.lagrange import Lagrange
    fiat_elem = Lagrange(ref_el, deg)

    my_elem = cg.to_fiat_elem()

    Q = create_quadrature(ref_el, 2*(deg+1))
    Qpts, _ = Q.get_points(), Q.get_weights()

    fiat_vals = fiat_elem.tabulate(0, Qpts)
    my_vals = my_elem.tabulate(0, Qpts)

    assert np.allclose(fiat_vals[(0,) * sd], my_vals[(0,) * sd])


@pytest.mark.parametrize("elem_gen,elem_code,deg", [(create_cg1, "CG", 1),
                                                    (create_dg1, "DG", 1),
                                                    (create_dg2, "DG", 2),
                                                    (create_cg2, "CG", 2)])
def test_entity_perms(elem_gen, elem_code, deg):
    cell = edge
    elem = elem_gen(cell)

    print(elem.to_fiat_elem())


# @pytest.mark.parametrize("elem_gen,elem_code,deg", [(create_cg1, "CG", 1),
#                                                     (create_dg1, "DG", 1),
#                                                     (create_dg2, "DG", 2),
#                                                     (create_cg2, "CG", 2)])
# def test_2d(elem_gen, elem_code, deg):
#     cell = edge
#     elem = elem_gen(cell)

#     mesh = UnitIntervalMesh(5)
#     V = FunctionSpace(mesh, elem_code, deg)
#     u = TrialFunction(V)
#     v = TestFunction(V)
#     f = Function(V)
#     x, = SpatialCoordinate(mesh)
#     f.interpolate(cos(x*pi*2))
#     a = (inner(grad(u), grad(v)) + inner(u, v)) * dx
#     L = inner(f, v) * dx
#     u1 = Function(V)
#     solve(a == L, u1)

#     V = FunctionSpace(mesh, elem.to_ufl_elem())
#     u = TrialFunction(V)
#     v = TestFunction(V)
#     f = Function(V)
#     x, = SpatialCoordinate(mesh)
#     f.interpolate((1+8*pi*pi)*cos(x*pi*2))
#     a = (inner(grad(u), grad(v)) + inner(u, v)) * dx
#     L = inner(f, v) * dx
#     u2 = Function(V)
#     solve(a == L, u2)

#     res = sqrt(assemble(dot(u1 - u1, u1 - u2) * dx))
#     assert np.allclose(res, 0)


# @pytest.mark.parametrize("elem_gen,elem_code,deg", [(create_cg1, "CG", 1),
#                                                     (create_dg1, "DG", 1),
#                                                     (construct_cg3, "CG", 3)])
# # ,
# #                                                     (lambda x:x, "CG", 4)
# def test_helmholtz(elem_gen, elem_code, deg):
#     cell = n_sided_polygon(3)
#     elem = elem_gen(cell)

#     mesh = UnitSquareMesh(20, 20)

#     V = FunctionSpace(mesh, elem_code, deg)
#     res1 = helmholtz_solve(mesh, V)

#     V = FunctionSpace(mesh, elem.to_ufl_elem())
#     res2 = helmholtz_solve(mesh, V)

#     res = sqrt(assemble(dot(res1 - res2, res1 - res2) * dx))
#     assert np.allclose(res, 0)


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
    return u


@pytest.mark.parametrize("cell", [vert, edge, tri])
def test_ufl_cell_conversion(cell):
    existing_cell = simplex(len(cell.vertices()))
    print(type(existing_cell))
    ufl_cell = cell.to_ufl()
    print(isinstance(ufl_cell, ufl.Cell))
    print(ufl_cell.cell_complex)
    print(ufl_cell.cellname())
