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
    Pk = PolynomialSpace(cell.dim(), cell.dim())
    dg = ElementTriple(cell, (Pk, CellL2, C0), DOFGenerator(xs, get_cyc_group(len(cell.vertices())), S1))
    return dg


def create_cg1(cell):
    deg = 1
    vert_dg = create_dg1(cell.vertices(get_class=True)[0])
    xs = [immerse(cell, vert_dg, TrH1)]

    Pk = PolynomialSpace(deg, deg)
    cg = ElementTriple(cell, (Pk, CellL2, C0), DOFGenerator(xs, get_cyc_group(len(cell.vertices())), S1))
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


@pytest.mark.parametrize("params", [(create_cg1, "CG", 1)])
# (create_dg1, "DG", 1)
def test_helmholtz(params):
    elem_gen, elem_code, deg = params
    cell = n_sided_polygon(3)
    elem = elem_gen(cell)

    mesh = UnitSquareMesh(10, 10)

    V = FunctionSpace(mesh, elem.to_ufl_elem())
    res = helmholtz_solve(mesh, V)

    V = FunctionSpace(mesh, elem_code, deg)
    res2 = helmholtz_solve(mesh, V)

    assert np.allclose(res, res2)


def helmholtz_solve(mesh, V):
    u = TrialFunction(V)
    v = TestFunction(V)
    f = Function(V)
    x, y = SpatialCoordinate(mesh)
    f.interpolate((1+8*pi*pi)*cos(x*pi*2)*cos(y*pi*2))
    a = (inner(grad(u), grad(v)) + inner(u, v)) * dx
    L = inner(f, v) * dx
    u = Function(V)
    solve(a == L, u, solver_parameters={'ksp_type': 'cg', 'pc_type': 'none'})
    f.interpolate(cos(x*pi*2)*cos(y*pi*2))
    res = sqrt(assemble(dot(u - f, u - f) * dx))
    return res


@pytest.mark.parametrize("cell", [vert, edge, tri])
def test_ufl_cell_conversion(cell):
    existing_cell = simplex(len(cell.vertices()))
    print(type(existing_cell))
    ufl_cell = cell.to_ufl()
    print(isinstance(ufl_cell, ufl.Cell))
    print(ufl_cell.cell_complex)
    print(ufl_cell.cellname())
