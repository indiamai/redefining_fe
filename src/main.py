# top level file for linking together all the packages
from firedrake import *
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


# set up of DG1
x = sp.Symbol("x")
v = sp.Function("v")(x)
# func1 = v(x)
xs = [lambda g: v.subs({"x": g(x)})]
dg1 = Triple(edges[0], ("L2"), E(xs, S2(), S1()))
ls = dg1.generate()
print(ls)

# def my_v(*x):
#     print("in my v")
#     print(x)
#     return -1
# from sympy import log, sin, cos, tan, Wild, Mul, Add
# from sympy.abc import y
# print(sp.sympify("k+x").evalf(subs={"k": 1}))
# f = sp.log(sp.sin(y)) + sp.tan(sp.sin(y**2))
# print(f.replace(sp.sin, sp.cos))
# print(f.replace(sp.sin, my_v))

# for l in ls:
#     print(l)
#     print(sp.srepr(l))
#     print(l.replace(v, my_v))
#     # print(l.replace(v(*), my_v(*)))
    
