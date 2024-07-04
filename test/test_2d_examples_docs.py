from redefining_fe import *
import pytest
import sympy as sp
import numpy as np


def test_dg_examples():
    # [test_dg0_definition]
    vert = Point(0)
    xs = [DOF(DeltaPairing(), PointKernel(()))]
    dg0 = ElementTriple(vert, (P0, CellL2, C0), DOFGenerator(xs, S1, S1))

    test_func = MyTestFunction(lambda: 3)

    for dof in dg0.generate():
        assert np.allclose(dof.eval(test_func), 3)

    edge = Point(1, [Point(0), Point(0)], vertex_num=2)

    xs = [DOF(DeltaPairing(), PointKernel((-1,)))]
    dg1 = ElementTriple(edge, (P1, CellL2, C0), DOFGenerator(xs, S2, S1))

    x = sp.Symbol("x")
    test_func = MyTestFunction(2*x, symbols=(x,))

    dof_vals = [-2, 2]

    for dof in dg1.generate():
        # Avoiding assuming order of generation
        assert any(np.isclose(val, dof.eval(test_func)) for val in dof_vals)

    tri = n_sided_polygon(3)
    xs = [DOF(DeltaPairing(), PointKernel((-1, -np.sqrt(3)/3)))]
    dg1 = ElementTriple(tri, (P1, CellL2, C0), DOFGenerator(xs, S3/S2, S1))

    dof_vals = [-11, 2, 9]

    x = sp.Symbol("x")
    y = sp.Symbol("y")
    test_func = MyTestFunction(10*x + 3*y/np.sqrt(3), symbols=(x, y))

    for dof in dg1.generate():
        assert any(np.isclose(val, dof.eval(test_func)) for val in dof_vals)


def test_cg_examples():
    # [test_cg_definition]
    tri = n_sided_polygon(3)
    edge = tri.edges(get_class=True)[0]
    vert = tri.vertices(get_class=True)[0]

    xs = [DOF(DeltaPairing(), PointKernel(()))]
    dg0 = ElementTriple(vert, (P0, CellL2, C0), DOFGenerator(xs, S1, S1))

    xs = [immerse(edge, dg0, CellH1)]
    cg1 = ElementTriple(edge, (P1, CellH1, C0),
                        DOFGenerator(xs, S2, S1))

    x = sp.Symbol("x")
    test_func = MyTestFunction(2*x, symbols=(x,))

    val_set = set([-2, 2])

    for dof in cg1.generate():
        # Avoiding assuming order of generation
        assert any(np.isclose(val, dof.eval(test_func)) for val in val_set)

    v_xs = [immerse(tri, dg0, CellH1)]
    v_dofs = DOFGenerator(v_xs, S3/S2, S1)

    xs = [DOF(DeltaPairing(), PointKernel((-1/3)))]
    dg0_int = ElementTriple(edge, (P1, CellH1, C0), DOFGenerator(xs, S2, S1))

    e_xs = [immerse(tri, dg0_int, CellH1)]
    e_dofs = DOFGenerator(e_xs, S3, S1)

    i_xs = [lambda g: DOF(DeltaPairing(), PointKernel(g((0, 0))))]
    i_dofs = DOFGenerator(i_xs, S1, S1)

    cg3 = ElementTriple(tri, (P3, CellH1, C0), [v_dofs, e_dofs, i_dofs])

    x = sp.Symbol("x")
    y = sp.Symbol("y")
    test_func = MyTestFunction(sp.Matrix([10*x, 3*y/np.sqrt(3)]), symbols=(x, y))

    dof_vals = np.array([[-10, -1], [0, 2], [10, -1],
                         [-10/3, 1], [-20/3, 0], [10/3, 1],
                         [20/3, 0], [-10/3, -1], [10/3, -1],
                         [0, 0]])

    for dof in cg3.generate():
        assert any([np.allclose(val, dof.eval(test_func).flatten()) for val in dof_vals])
