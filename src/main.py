# top level file for linking together all the packages
from firedrake import *
import numpy as np
from groups.new_groups import r, rot, S1, S2, S3, D4
from cell_complex.cells import Point, Edge
from FIAT.functional import PointEvaluation
from FIAT.reference_element import Point as fiatPoint, UFCInterval, UFCTriangle
from triples import Triple, E, immerse
from spaces.continuous_spaces import H1, L2



vertices = []
for i in range(3):
    vertices.append(Point(0))
edges = []
edges.append(
    Point(1, [Edge(vertices[0], lambda x: -1),
              Edge(vertices[1], lambda x: 1)]))
edges.append(
    Point(1, [Edge(vertices[0], lambda x: 1),
              Edge(vertices[2], lambda x: -1)]))
edges.append(
    Point(1, [Edge(vertices[1], lambda x: -1),
              Edge(vertices[2], lambda x: 1)]))

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

# e2 = edges[0].orient(r)
# print("original")
# print(edges[0].G)
# print(edges[0].graph().edges)
# print(edges[0].graph()[3][1]["edge_class"])
# edge_class1 = edges[0].cell_attachment_route(0)


# print("copied")
# print(e2.G)
# print(e2.G[3][1]["edge_class"])
# edge_class2 = e2.cell_attachment_route(0)
# edges[0].plot()
# e2.plot()


# dg0 on point
print("DG0 on point")
ref_elem = fiatPoint()
xs = [lambda g: PointEvaluation(ref_elem, g(()))]
dg0 = Triple(vertices[0], ("P0", L2(), "C0"), E(xs, S1, S1))
ls = dg0.generate()
for dof in ls:
    print(dof.tostr())

# cg1 on interval
print("CG1 on interval")
xs = [lambda g: immerse(g, edges[0].cell_attachment(0), dg0)]
cg1 = Triple(edges[0], ("P1", H1(), "C0"), E(xs, S2, S1))
ls = cg1.generate()
for dof in ls:
    print(dof.tostr())

# # dg1 on interval
print("DG1 on interval")
ref_interval = UFCInterval()
xs = [lambda g: PointEvaluation(ref_interval, g((-1,)))]
dg1 = Triple(edges[0], ("P1", "L2", "C0"), E(xs, S2, S1))
ls = dg1.generate()
for dof in ls:
    print(dof.tostr())

print((S3/S2).members)
# dg1 on triangle
print("DG1 on triangle")
ref_triangle = UFCTriangle()
xs = [lambda g: PointEvaluation(ref_triangle, g((-1, -np.sqrt(3)/3)))]
dg1 = Triple(a4, ("P2", "L2", "C0"), E(xs, S3/S2, S1))
ls = dg1.generate()
for dof in ls:
    print(dof.tostr())

print((D4/S2).members)
