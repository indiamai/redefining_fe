# top level file for linking together all the packages
from firedrake import *
import numpy as np
import sympy as sp
from groups.groups import r, rot, S1, S2, S3, D4, C4, S4
from cell_complex.cells import Point, Edge, n_sided_polygon
from dof_lang.dof import DeltaPairing, DOF, L2InnerProd, MyTestFunction, PointKernel, PolynomialKernel
from triples import ElementTriple, DOFGenerator, immerse
from spaces.element_sobolev_spaces import CellH1, CellL2, CellHDiv, CellHCurl, CellH2, CellH3
from spaces.polynomial_spaces import P0, P1, P2, P3, Q2, VectorPolynomialSpace
from spaces.interpolation_spaces import C0, L2, H1, HDiv
import matplotlib.pyplot as plt
from ufl.sobolevspace import L2 as uflL2, H1 as uflH1, HDiv as uflHDiv
# import json
# import jsonpickle

# from scratch import n_sided_polygon

vertices = []
for i in range(2):
    vertices.append(Point(0))
edges = []
edges.append(
    Point(1, [Edge(vertices[0], lambda: (-1,)),
              Edge(vertices[1], lambda: (1,))]))
# edges.append(
#     Point(1, [Edge(vertices[0], lambda: (1,)),
#               Edge(vertices[2], lambda: (-1,))]))
# edges.append(
#     Point(1, [Edge(vertices[1], lambda: (-1,)),
#               Edge(vertices[2], lambda: (1,))]))

# a4 = Point(2, [Edge(edges[0], lambda x: [x, -np.sqrt(3) / 3]),
#                Edge(edges[1], lambda x: [(- x - 1) / 2,
#                                          np.sqrt(3) * (3 * -x + 1) / 6]),
#                Edge(edges[2], lambda x: [(1 - x) / 2,
#                                          np.sqrt(3) * (3 * x + 1) / 6])])
# print(a4.vertices())
# print(a4.basis_vectors())
# print(a4.basis_vectors(return_coords=True))
# print(a4.basis_vectors(entity=edges[0], return_coords=True))
# print(a4.basis_vectors(entity=edges[1], return_coords=True))
# print(a4.basis_vectors(entity=edges[2], return_coords=True))
# print(a4.basis_vectors(entity=edges[0]))
# print(a4.basis_vectors(entity=edges[1]))
# print(a4.basis_vectors(entity=edges[2]))
# a4rot = a4.orient(rot)
# # a4.hasse_diagram()
# # a4.plot()
# # a4rot.plot()
# # edges[0].plot()
# # e2 = edges[0].orient(r)
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

# intervalH1 = CellH1(edges[0])
# intervalHDiv = CellHDiv(edges[0])
# intervalHCurl = CellHCurl(edges[0])
# pointL2 = CellL2(vertices[0])
# intervalL2 = CellL2(edges[0])
# triangleL2 = CellL2(a4)
# triHCurl = CellHCurl(a4)
# triHDiv = CellHDiv(a4)
# triangleH1 = CellH1(a4)
# triangleH2 = CellH2(a4)
# print(intervalH1)
# print(intervalHDiv)
# # this comparision should also check that the cells are the same (subspaces?)
# print(intervalH1 < intervalHDiv)


# # def test_p1_v(x):
# #     return 2*x + 3

x = sp.Symbol("x")

RTspace1 = P0 + x*P0
print(RTspace1.weights)
print(RTspace1.spaces)
print(RTspace1)
# edges[0].get_topology()
polyset = P0.to_ON_polynomial_set(edges[0])
quit()

tri = n_sided_polygon(3)
vert = tri.d_entities(0, get_class=True)[0]
edge = tri.d_entities(1, get_class=True)[0]


test_func = MyTestFunction(lambda x: 2*x + 3)
# # test_func2 = MyTestFunction(lambda x, y: (10*x, y))
# # print(test_func)
# # # dg0 on point
print("DG0 on point")
xs = [DOF(DeltaPairing(), PointKernel(()))]
dg0 = ElementTriple(vert, (P0, CellL2, C0),
                    DOFGenerator(xs, S1, S1))
ls = dg0.generate()
# print("num dofs ", dg0.num_dofs())
# for dof in ls:
#     print(dof)
# print(dg0.__dict__)
# json.dumps(dg0.__dict__)
# json_string = jsonpickle.encode(dg0)
# print(json_string)

# # cg1 on interval
print("CG1 on interval")
xs = [immerse(edge, dg0, CellH1)]
cg1 = ElementTriple(edge, (P1, CellH1, C0),
                    DOFGenerator(xs, S2, S1))
ls = cg1.generate()
# print("num dofs ", cg1.num_dofs())
# for dof in ls:
#     print(dof)
#     print(dof.eval(test_func))
# cg1.plot()

# # # # # dg1 on interval
print("DG1 on interval")
xs = [DOF(DeltaPairing(), PointKernel((-1,)))]
dg1 = ElementTriple(edge, (P1, CellL2, C0),
                    DOFGenerator(xs, S2, S1))
ls = dg1.generate()
# print("num dofs ", dg1.num_dofs())
# for dof in ls:
#     print(dof)
#     print(dof.eval(test_func))
# dg1.plot()
# # dg1 on triangle
print("DG1 on triangle")
xs = [DOF(DeltaPairing(), PointKernel((-1, -np.sqrt(3)/3)))]
dg1 = ElementTriple(tri, (P1, CellL2, C0),
                    DOFGenerator(xs, S3/S2, S1))
ls = dg1.generate()
# print("num dofs ", dg1.num_dofs())
# for dof in ls:
#     print(dof)

# print("DG0 on interval")
xs = [DOF(DeltaPairing(), PointKernel((0,)))]
dg0_int = ElementTriple(edge, (P0, CellL2, C0),
                        DOFGenerator(xs, S1, S1))
ls = dg0_int.generate()
# print("num dofs ", dg0_int.num_dofs())
# for dof in ls:
#     print(dof)
# dg0_int.plot()

# # # cg3 on triangle
print("CG3")
v_xs = [immerse(tri, dg0, CellH1)]
v_dofs = DOFGenerator(v_xs, S3/S2, S1)

e_xs = [immerse(tri, dg0_int, CellH1)]
e_dofs = DOFGenerator(e_xs, S3, S1)

i_xs = [lambda g: DOF(DeltaPairing(), PointKernel(g((0, 0))))]
i_dofs = DOFGenerator(i_xs, S1, S1)

cg3 = ElementTriple(tri, (P3, CellH1, C0),
                    [v_dofs, e_dofs, i_dofs])

phi_0 = MyTestFunction(lambda x, y: (x, y))
ls = cg3.generate()
# print("num dofs ", cg3.num_dofs())
# for dof in ls:
#     print(dof)
#     print(dof.eval(phi_0))
# cg3.plot()

print("Integral Moment")
xs = [DOF(L2InnerProd(), PointKernel((1,)))]
dofs = DOFGenerator(xs, S1, S2)

int_ned = ElementTriple(edge, (P1, CellHCurl, C0), dofs)
ls = int_ned.generate()
# for dof in ls:
#     print(dof)

phi_2 = MyTestFunction(lambda x, y: (1/3 - (np.sqrt(3)/6)*y,
                                     (np.sqrt(3)/6)*x))
phi_0 = MyTestFunction(lambda x, y: (-1/6 - (np.sqrt(3)/6)*y,
                                     (-np.sqrt(3)/6) + (np.sqrt(3)/6)*x))
phi_1 = MyTestFunction(lambda x, y: (-1/6 - (np.sqrt(3)/6)*y,
                                     (np.sqrt(3)/6) + (np.sqrt(3)/6)*x))

print("Nedelec")
xs = [immerse(tri, int_ned, CellHCurl)]
tri_dofs = DOFGenerator(xs, S3, S3)
vecP3 = VectorPolynomialSpace(P3, P3)
ned = ElementTriple(tri, (vecP3, CellHCurl, C0), [tri_dofs])
# ned.plot()
ls = ned.generate()
# for dof in ls:
#     print(dof)
#     print("phi_0 ", dof.eval(phi_0))
#     print("phi_1 ", dof.eval(phi_1))
#     print("phi_2 ", dof.eval(phi_2))

print("Edge of RT")
xs = [DOF(L2InnerProd(), PointKernel((1,)))]
dofs = DOFGenerator(xs, S1, S2)
int_rt = ElementTriple(edge, (P1, CellHDiv, C0), dofs)
ls = int_rt.generate()
# for dof in ls:
#     print(dof)

print("RT")

phi_2 = MyTestFunction(lambda x, y: ((np.sqrt(3)/6)*x,
                                     -1/3 + (np.sqrt(3)/6)*y))
phi_0 = MyTestFunction(lambda x, y: ((-np.sqrt(3)/6) + (np.sqrt(3)/6)*x,
                                     1/6 + (np.sqrt(3)/6)*y))
phi_1 = MyTestFunction(lambda x, y: ((np.sqrt(3)/6) + (np.sqrt(3)/6)*x,
                                     1/6 + (np.sqrt(3)/6)*y))

xs = [immerse(tri, int_rt, CellHDiv)]
tri_dofs = DOFGenerator(xs, S3, S3)
vecP3 = VectorPolynomialSpace(P3, P3)
rt = ElementTriple(tri, (vecP3, CellHDiv, C0), [tri_dofs])
# ls = rt.generate()
# for dof in ls:
#     print(dof)
#     print("phi_0 ", dof.eval(phi_0))
#     print("phi_1 ", dof.eval(phi_1))
#     print("phi_2 ", dof.eval(phi_2))
# rt.plot()

print("Hermite")
v_xs = [immerse(tri, dg0, CellH1)]
v_dofs = DOFGenerator(v_xs, S3/S2, S1)

v_derv_xs = [immerse(tri, dg0, CellH2)]
v_derv_dofs = DOFGenerator(v_derv_xs, S3/S2, S1)

v_derv2_xs = [immerse(tri, dg0, CellH3)]
v_derv2_dofs = DOFGenerator(v_derv2_xs, S3/S2, S1)

i_xs = [DOF(DeltaPairing(), PointKernel((0, 0)))]
i_dofs = DOFGenerator(i_xs, S1, S1)

her = ElementTriple(tri, (P3, CellH2, C0),
                    [v_dofs, v_derv_dofs, v_derv2_dofs, i_dofs])

phi_0 = MyTestFunction(lambda x, y: x**2 + 3*y**3 + 4*x*y)
ls = her.generate()
# print("num dofs ", her.num_dofs())
# for dof in ls:
#     print(dof)
# #     # print("dof eval", dof.eval(phi_0))
    
# her.plot()


print("CG1 on Square")
square = n_sided_polygon(4)
# tri2.plot()
# print(square.vertices(return_coords=True))
# print(square.group.size())
vert = square.d_entities(0, get_class=True)[0]
edge = square.d_entities(1, get_class=True)[0]

xs = [DOF(DeltaPairing(), PointKernel(()))]
dg0 = ElementTriple(vert, (P0, CellL2, C0),
                    DOFGenerator(xs, S1, S1))

xs = [DOF(DeltaPairing(), PointKernel((0,)))]
dg0_int = ElementTriple(edge, (P0, CellL2, C0),
                        DOFGenerator(xs, S1, S1))

print("CG3")
v_xs = [immerse(square, dg0, CellH1)]
v_dofs = DOFGenerator(v_xs, C4, S1)

e_xs = [immerse(square, dg0_int, CellH1)]
e_dofs = DOFGenerator(e_xs, D4, S1)

i_xs = [lambda g: DOF(DeltaPairing(), PointKernel(g((0, 0))))]
i_dofs = DOFGenerator(i_xs, S1, S1)

cg3 = ElementTriple(square, (P3, CellH1, C0),
                    [v_dofs, e_dofs, i_dofs])

phi_0 = MyTestFunction(lambda x, y: (x, y))
ls = cg3.generate()
# print("num dofs ", cg3.num_dofs())
# for dof in ls:
#     print(dof)
#     print(dof.eval(phi_0))
# cg3.plot()

print("Edge of RT 2nd order")
xs = [DOF(L2InnerProd(), PolynomialKernel(lambda x: (1/2)*(1 + x)))]
dofs = DOFGenerator(xs, S2, S2)
int_rt2 = ElementTriple(edge, (P1, CellHDiv, C0), dofs)
# ls = int_rt2.generate()
# for l in ls:
#     print(l.eval(MyTestFunction(lambda x: x)))
    # print(l.eval(MyTestFunction(lambda x: 0)))


xs = [immerse(tri, int_rt2, CellHDiv)]
tri_dofs = DOFGenerator(xs, S3, S3)

i_xs = [lambda g: DOF(L2InnerProd(), PointKernel(g((1, 0)))),
        lambda g: DOF(L2InnerProd(), PointKernel(g((0, 1))))]
i_dofs = DOFGenerator(i_xs, S1, S1)

vecP3 = VectorPolynomialSpace(P3, P3)
rt2 = ElementTriple(tri, (vecP3, CellHDiv, C0), [tri_dofs, i_dofs])
ls = rt2.generate()
# for l in ls:
#     print(l)

