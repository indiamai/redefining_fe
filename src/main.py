# top level file for linking together all the packages
# from firedrake import *
import numpy as np
import sympy as sp
from groups.groups import r, rot, S1, S2, S3
from cell_complex.cells import Point, Edge
from triples import Triple, E


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
# a4.plot()
# print(a4.vertices())
# print(a4.basis_vectors(return_coords=True))
a4.orient(rot)
# print(a4.basis_vectors(return_coords=True))
# a4.hasse_diagram()
# a4.plot()


# set up of CG1
x = sp.Symbol("x")
v = sp.Function("v")(x)
x = [lambda g: v.subs({"x": g(x)})]
dg1 = Triple(edges[0], ("L2"), E(x, S2(), S1()))
print(dg1.generate())
