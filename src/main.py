# top level file for linking together all the packages
from firedrake import *
import numpy as np
from groups.new_groups import r, rot, S1, S2, S3, D4, C4
from cell_complex.cells import Point, Edge
from dof_lang.dof import construct_point_eval, construct_tangent_dof, MyTestFunction
from triples import ElementTriple, DOFGenerator, immerse
from spaces.element_sobolev_spaces import CellH1, CellL2, CellHDiv, CellHCurl
from spaces.polynomial_spaces import P0, P1, P2, P3, Q2


vertices = []
for i in range(3):
    vertices.append(Point(0))
edges = []
edges.append(
    Point(1, [Edge(vertices[0], lambda: (-1,)),
              Edge(vertices[1], lambda: (1,))]))
edges.append(
    Point(1, [Edge(vertices[0], lambda: (1,)),
              Edge(vertices[2], lambda: (-1,))]))
edges.append(
    Point(1, [Edge(vertices[1], lambda: (-1,)),
              Edge(vertices[2], lambda: (1,))]))

a4 = Point(2, [Edge(edges[0], lambda x: [x, -np.sqrt(3) / 3]),
               Edge(edges[1], lambda x: [(- x - 1) / 2,
                                         np.sqrt(3) * (3 * -x + 1) / 6]),
               Edge(edges[2], lambda x: [(1 - x) / 2,
                                         np.sqrt(3) * (3 * x + 1) / 6])])
# print(a4.vertices())
# print(a4.basis_vectors())
# print(a4.basis_vectors(return_coords=True))
# a4rot = a4.orient(rot)
# a4.hasse_diagram()
# a4.plot()
# a4rot.plot()
# edges[0].plot()
# e2 = edges[0].orient(r)
# print("original")
# print(edges[0].G)
# print(edges[0].graph().edges)
# print(edges[0].graph()[3][1]["edge_class"])
# edge_class1 = edges[0].cell_attachment(0)


# print("copied")
# print(e2.G)
# print(e2.G[3][1]["edge_class"])
# edge_class2 = e2.cell_attachment_route(0)
# edges[0].plot()
# e2.plot()

intervalH1 = CellH1(edges[0])
intervalHDiv = CellHDiv(edges[0])
intervalHCurl = CellHCurl(edges[0])
pointL2 = CellL2(vertices[0])
intervalL2 = CellL2(edges[0])
triangleL2 = CellL2(a4)
triHCurl = CellHCurl(a4)
triangleH1 = CellH1(a4)
print(intervalH1)
print(intervalHDiv)
# this comparision should also check that the cells are the same (subspaces?)
print(intervalH1 < intervalHDiv)


# def test_p1_v(x):
#     return 2*x + 3


test_func = MyTestFunction(lambda x: 2*x + 3)
test_func2 = MyTestFunction(lambda x, y: (10*x, y))
print(test_func)
# dg0 on point
print("DG0 on point")
xs = [lambda g: construct_point_eval(g(()), vertices[0], pointL2)]
dg0 = ElementTriple(vertices[0], (P0, pointL2, "C0"),
                    DOFGenerator(xs, S1, S1))
ls = dg0.generate()
print("num dofs ", dg0.num_dofs())
for dof in ls:
    print(dof)

# cg1 on interval
print("CG1 on interval")
xs = [lambda g: immerse(g, edges[0], dg0, intervalH1)]
cg1 = ElementTriple(edges[0], (P1, intervalH1, "C0"),
                    DOFGenerator(xs, S2, S1))
ls = cg1.generate()
print("num dofs ", cg1.num_dofs())
for dof in ls:
    print(dof)
    print(dof(test_func))

# # dg1 on interval
print("DG1 on interval")
xs = [lambda g: construct_point_eval(g((-1,)), edges[0], intervalL2)]
dg1 = ElementTriple(edges[0], (P1, intervalL2, "C0"),
                    DOFGenerator(xs, S2, S1))
ls = dg1.generate()
print("num dofs ", dg1.num_dofs())
for dof in ls:
    print(dof)
    print(dof(test_func))

# dg1 on triangle
print("DG1 on triangle")
xs = [lambda g: construct_point_eval(g((-1, -np.sqrt(3)/3)), a4, triangleL2)]
dg1 = ElementTriple(a4, (P1, triangleL2, "C0"),
                    DOFGenerator(xs, S3/S2, S1))
ls = dg1.generate()
print("num dofs ", dg1.num_dofs())
for dof in ls:
    print(dof)

print("DG0 on interval")
xs = [lambda g: construct_point_eval(g((0,)), edges[0], intervalL2)]
dg0_int = ElementTriple(edges[0], (P0, intervalL2, "C0"),
                        DOFGenerator(xs, S1, S1))
ls = dg0_int.generate()
print("num dofs ", dg0_int.num_dofs())
for dof in ls:
    print(dof)
# dg0_int.plot()

# cg3 on triangle
print("CG3")
v_xs = [lambda g: immerse(g, a4, dg0, triangleH1)]
v_dofs = DOFGenerator(v_xs, S3/S2, S1)

e_xs = [lambda g: immerse(g, a4, dg0_int, triangleH1)]
e_dofs = DOFGenerator(e_xs, S3/S2, S1)

i_xs = [lambda g: construct_point_eval(g((0, 0)), a4, triangleL2)]
i_dofs = DOFGenerator(i_xs, S1, S1)

cg3 = ElementTriple(a4, (P3, triangleH1, "C0"),
                    [v_dofs, e_dofs, i_dofs])

ls = cg3.generate()
print("num dofs ", cg3.num_dofs())
for dof in ls:
    print(dof)
    # print(dof(test_func2))
cg3.plot()


# print("Integral Moment")
# xs = [lambda g: construct_tangent_dof(g(edges[0]), intervalHCurl)]
# dofs = DOFGenerator(xs, S1, S2)

# int_ned = ElementTriple(edges[0], (P1, intervalHCurl, "C0"), dofs)
# ls = int_ned.generate()
# for dof in ls:
#     print(dof)

# print("Rotation of integral moments")
# xs = [lambda g: construct_tangent_dof(g(a4).edges(get_class=True)[0], triHCurl)]
# tri_dofs = DOFGenerator(xs, S3/S2, S3)

# ned_dg = ElementTriple(a4, (P3, triHCurl, "C0"),
#                       [tri_dofs])
# ls = ned_dg.generate()
# for dof in ls:
#     print(dof)

# xs = [lambda g: immerse(g, a4, int_ned, triHCurl)]
# tri_dofs = DOFGenerator(xs, S3/S2, S3)

# ned = ElementTriple(a4, (P3, triHCurl, "C0"),
#                     [tri_dofs])
# ls = ned.generate()
# for dof in ls:
#     print(dof)