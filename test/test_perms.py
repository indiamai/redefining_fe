from redefining_fe import *
from test_convert_to_fiat import create_cg1, create_dg1, create_cg2
from test_2d_examples_docs import construct_cg3
import pytest
import sympy as sp
import numpy as np

vert = Point(0)
edge = Point(1, [Point(0), Point(0)], vertex_num=2)
tri = n_sided_polygon(3)


@pytest.mark.parametrize("cell", [edge])
def test_cg_perms(cell):
    cg1 = create_cg1(cell)
    cg1.make_dof_perms()


@pytest.mark.parametrize("cell", [edge])
def test_cg2_perms(cell):
    cg2 = create_cg2(cell)

    cg2.make_dof_perms()
    for g in cg2.generate():
        print(g)
        print(g.generation)


@pytest.mark.parametrize("cell", [tri])
def test_cg3_perms(cell):
    cg3 = construct_cg3(cell)
    cg3.make_dof_perms()


@pytest.mark.parametrize("cell", [vert, edge])
def test_dg_perms(cell):
    cg1 = create_dg1(cell)
    cg1.make_dof_perms()


def test_nd_perms():
    xs = [DOF(L2InnerProd(), PointKernel(edge.basis_vectors()[0]))]
    dofs = DOFGenerator(xs, S1, S2)
    int_ned = ElementTriple(edge, (P1, CellHCurl, C0), dofs)
    int_ned.make_dof_perms()


@pytest.mark.parametrize("cell", [edge, tri])
def test_basic_perms(cell):
    cell_group = cell.group
    bvs = cell.basis_vectors()
    M = np.array(bvs).T
    for g in cell_group.members():
        trans_bvs = np.array([g(bvs[i]) for i in range(len(bvs))]).T
        print(np.linalg.solve(M, trans_bvs))


def test_square():
    square = n_sided_polygon(4)
    edge = square.d_entities(1, get_class=True)[0]

    xs = [DOF(L2InnerProd(), PointKernel(edge.basis_vectors()[0]))]
    dg0_int = ElementTriple(edge, (P0, CellL2, C0),
                            DOFGenerator(xs, S1, S2))

    e_xs = [immerse(square, dg0_int, TrH1)]
    e_dofs = DOFGenerator(e_xs, sq_edges, S1)

    i_xs = [lambda g: DOF(DeltaPairing(), PointKernel(g((0, 0))))]
    i_dofs = DOFGenerator(i_xs, S1, S1)

    cg3 = ElementTriple(square, (P3, CellH1, C0),
                        [e_dofs, i_dofs])
    for dof in cg3.generate():
        print(dof)
    cg3.make_dof_perms()