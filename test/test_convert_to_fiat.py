import pytest
import numpy as np
from redefining_fe import *
from FIAT.quadrature_schemes import create_quadrature
from firedrake import *
from ufl.cell import simplex
# from firedrake import functionspaceimpl as impl
# import finat
# from FInAT.fiat_elements import FiatElement

vert = Point(0)
edge = Point(1, [Point(0), Point(0)], vertex_num=2)
tri = n_sided_polygon(3)


# @pytest.mark.parametrize("cell", [vert, edge, tri])
def create_dg1(cell):
    xs = [DOF(DeltaPairing(), PointKernel(cell.vertices(return_coords=True)[0]))]
    Pk = PolynomialSpace(cell.dim, cell.dim)
    dg = ElementTriple(cell, (Pk, CellL2, C0), DOFGenerator(xs, get_cyc_group(len(cell.vertices())), S1))
    return dg


@pytest.mark.parametrize("cell", [tri])
def test_create_cg1(cell):
    deg = 1
    vert_dg = create_dg1(cell.vertices(get_class=True)[0])
    xs = [immerse(cell, vert_dg, TrH1)]

    Pk = PolynomialSpace(deg, deg)
    cg = ElementTriple(cell, (Pk, CellL2, C0), DOFGenerator(xs, get_cyc_group(len(cell.vertices())), S1))

    from FIAT.lagrange import Lagrange
    ref_el = cell.to_fiat()
    sd = ref_el.get_spatial_dimension()

    fiat_elem = Lagrange(ref_el, deg)
    my_elem = cg.to_fiat_elem()

    Q = create_quadrature(ref_el, 2*(deg+1))
    Qpts, _ = Q.get_points(), Q.get_weights()

    fiat_vals = fiat_elem.tabulate(0, Qpts)
    my_vals = my_elem.tabulate(0, Qpts)

    assert np.allclose(fiat_vals[(0,) * sd], my_vals[(0,) * sd])
    mesh = UnitSquareMesh(2, 2)

    V = FunctionSpace(mesh, "CG", 1)

    u = TrialFunction(V)
    v = TestFunction(V)
    a = (inner(grad(u), grad(v)) + inner(u, v)) * dx

    a = assemble(a)

    ufl_elem = cg.to_ufl_elem()
    V_n = FunctionSpace(mesh, ufl_elem)
    u_n = TrialFunction(V_n)
    v_n = TestFunction(V_n)
    # a = (inner(grad(u_n), grad(v_n)) + inner(u_n, v_n)) * dx
    a = u_n * v_n * dx

    # a_n = assemble(a)


@pytest.mark.parametrize("cell", [vert, edge, tri])
def test_ufl_cell_conversion(cell):
    existing_cell = simplex(len(cell.vertices()))
    print(type(existing_cell))
    ufl_cell = cell.to_ufl()
    print(isinstance(ufl_cell, ufl.Cell))
    print(ufl_cell.fe_cell)
    print(ufl_cell.cellname())
