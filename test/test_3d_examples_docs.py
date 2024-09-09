from redefining_fe import *
import sympy as sp
import numpy as np


def test_dg1():
    tetra = make_tetrahedron()

    xs = [DOF(DeltaPairing(), PointKernel(tuple(tetra.vertices(return_coords=True)[0])))]
    dg1 = ElementTriple(tetra, (P1, CellL2, "C0"),
                        DOFGenerator(xs, Z4, S1))

    x = sp.Symbol("x")
    y = sp.Symbol("y")
    z = sp.Symbol("z")
    test_func = MyTestFunction(x + y + z, symbols=(x, y, z))

    dof_vals = [x+y+z for (x, y, z) in tetra.vertices(return_coords=True)]

    for dof in dg1.generate():
        assert any(np.isclose(val, dof.eval(test_func)) for val in dof_vals)


def test_tet_cg3():
    tetra = make_tetrahedron()
    vert = tetra.vertices(get_class=True)[0]
    edge = tetra.edges(get_class=True)[0]
    face = tetra.d_entities(2, get_class=True)[0]

    xs = [DOF(DeltaPairing(), PointKernel(()))]
    dg0 = ElementTriple(vert, (P0, CellL2, "C0"),
                        DOFGenerator(xs, S1, S1))

    xs = [DOF(DeltaPairing(), PointKernel((-1/3,)))]
    dg1_int = ElementTriple(edge, (P0, CellL2, "C0"),
                            DOFGenerator(xs, S2, S1))

    xs = [DOF(DeltaPairing(), PointKernel((0, 0)))]
    dg0_face = ElementTriple(face, (P0, CellL2, "C0"),
                             DOFGenerator(xs, S1, S1))

    v_xs = [immerse(tetra, dg0, TrH1)]
    cgverts = DOFGenerator(v_xs, Z4, S1)

    e_xs = [immerse(tetra, dg1_int, TrH1)]
    cgedges = DOFGenerator(e_xs, A4, S1)

    f_xs = [immerse(tetra, dg0_face, TrH1)]
    cgfaces = DOFGenerator(f_xs, S4, S1)

    cg3 = ElementTriple(tetra, (P1, CellH1, "C0"),
                        [cgverts, cgedges, cgfaces])

    x = sp.Symbol("x")
    y = sp.Symbol("y")
    z = sp.Symbol("z")
    test_func = MyTestFunction(sp.Matrix([10*x, 3*y/np.sqrt(3), z*4]), symbols=(x, y, z))

    for dof in cg3.generate():
        dof.eval(test_func)


def test_tet_rt():
    tetra = make_tetrahedron()
    face = tetra.d_entities(2, get_class=True)[0]

    xs = [DOF(L2InnerProd(), PolynomialKernel(lambda x: 1))]
    dofs = DOFGenerator(xs, S1, S2)
    face_vec = ElementTriple(face, (P1, CellHDiv, "C0"), dofs)
    ls = face_vec.generate()

    im_xs = [immerse(tetra, face_vec, TrHDiv)]
    face = DOFGenerator(im_xs, S4, S4)

    rt1 = ElementTriple(tetra, (P1, CellHDiv, "C0"),
                        [face])
    ls = rt1.generate()
    # TODO make this a proper test
    for dof in ls:
        print(dof)


def test_tet_ned():
    tetra = make_tetrahedron()
    edge = tetra.edges(get_class=True)[0]

    xs = [DOF(L2InnerProd(), PolynomialKernel(lambda x: 1))]
    dofs = DOFGenerator(xs, S1, S2)
    int_ned = ElementTriple(edge, (P1, CellHCurl, "C0"), dofs)
    ls = int_ned.generate()

    im_xs = [immerse(tetra, int_ned, TrHCurl)]
    edge = DOFGenerator(im_xs, A4, Z4)

    ned = ElementTriple(tetra, (P1, CellHCurl, "C0"),
                        [edge])
    ls = ned.generate()
    # TODO make this a proper test
    for dof in ls:
        print(dof)
