import pytest
import numpy as np
from redefining_fe import *
from FIAT.quadrature_schemes import create_quadrature
from firedrake import *
from ufl.cell import simplex
from test_2d_examples_docs import construct_cg3
# from firedrake import functionspaceimpl as impl
# import finat
# from FInAT.fiat_elements import FiatElement

vert = Point(0)
edge = Point(1, [Point(0), Point(0)], vertex_num=2)
tri = n_sided_polygon(3)


# @pytest.mark.parametrize("cell", [vert, edge, tri])
def create_dg1(cell):
    xs = [DOF(DeltaPairing(), PointKernel(cell.vertices(return_coords=True)[0]))]
    Pk = PolynomialSpace(cell.dim(), cell.dim())
    dg = ElementTriple(cell, (Pk, CellL2, C0), DOFGenerator(xs, get_cyc_group(len(cell.vertices())), S1))
    return dg


# @pytest.mark.parametrize("cell", [tri])
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
    print(type(mesh))


    ufl_elem = cg.to_ufl_elem()
    V_temp = VectorFunctionSpace(mesh, ufl_elem)
    new_mesh = Mesh(Function(V_temp))
    print(new_mesh.ufl_cell())
    V_n = FunctionSpace(new_mesh, ufl_elem)
    
    u_n = TrialFunction(V_n)
    v_n = TestFunction(V_n)
    a = u_n * v_n * dx

    a_n = assemble(a)

def test_helmholtz():
    cell = n_sided_polygon(3)
    deg = 1
    vert_dg = create_dg1(cell.vertices(get_class=True)[0])
    xs = [immerse(cell, vert_dg, TrH1)]

    Pk = PolynomialSpace(deg, deg)
    cg = ElementTriple(cell, (Pk, CellL2, C0), DOFGenerator(xs, get_cyc_group(len(cell.vertices())), S1))
    dg = create_dg1(cell)
    # cg3 = construct_cg3()
    mesh = UnitSquareMesh(10, 10)
    print(type(mesh.ufl_cell()))
    V = FunctionSpace(mesh, dg.to_ufl_elem())
    u = TrialFunction(V)
    v = TestFunction(V)
    f = Function(V)
    x, y = SpatialCoordinate(mesh)
    f.interpolate((1+8*pi*pi)*cos(x*pi*2)*cos(y*pi*2))
    print(type(u.function_space().mesh().ufl_cell()))
    a = (inner(grad(u), grad(v)) + inner(u, v)) * dx
    L = inner(f, v) * dx

    u = Function(V)

    solve(a == L, u, solver_parameters={'ksp_type': 'cg', 'pc_type': 'none'})

    f.interpolate(cos(x*pi*2)*cos(y*pi*2))
    print(sqrt(assemble(dot(u - f, u - f) * dx)))

@pytest.mark.parametrize("cell", [vert, edge, tri])
def test_ufl_cell_conversion(cell):
    existing_cell = simplex(len(cell.vertices()))
    print(type(existing_cell))
    ufl_cell = cell.to_ufl()
    print(isinstance(ufl_cell, ufl.Cell))
    print(ufl_cell.fe_cell)
    print(ufl_cell.cellname())
