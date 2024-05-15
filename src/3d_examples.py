# Examples of Elements in 3 dimensions
from firedrake import *
import numpy as np
from sympy.combinatorics import PermutationGroup, Permutation
from sympy.combinatorics.named_groups import SymmetricGroup, DihedralGroup, CyclicGroup, AlternatingGroup
from groups.new_groups import r, r_y, rot, S1, S2, S3, D4, C3, C4, S4, GroupRepresentation
from cell_complex.cells import Point, Edge
from dof_lang.dof import DeltaPairing, DOF, L2InnerProd, MyTestFunction, PointKernel
from triples import ElementTriple, DOFGenerator, immerse
from spaces.element_sobolev_spaces import CellH1, CellL2, CellHDiv, CellHCurl, CellH2, CellH3
from spaces.polynomial_spaces import P0, P1, P2, P3, Q2, VectorPolynomialSpace
import matplotlib.pyplot as plt


vertices = []
for i in range(4):
    vertices.append(Point(0))

edges = []
edges.append(
        Point(1, [Edge(vertices[0], lambda: (-1,)),
                  Edge(vertices[1], lambda: (1,))]))
edges.append(
        Point(1, [Edge(vertices[1], lambda: (-1,)),
                  Edge(vertices[2], lambda: (1,))]))
edges.append(
        Point(1, [Edge(vertices[2], lambda: (-1,)),
                  Edge(vertices[0], lambda: (1,))]))
edges.append(
        Point(1, [Edge(vertices[3], lambda: (-1,)),
                  Edge(vertices[0], lambda: (1,))]))
edges.append(
        Point(1, [Edge(vertices[1], lambda: (-1,)),
                  Edge(vertices[3], lambda: (1,))]))
edges.append(
        Point(1, [Edge(vertices[2], lambda: (-1,)),
                  Edge(vertices[3], lambda: (1,))]))

face1 = Point(2, [Edge(edges[3], lambda x: [x, -np.sqrt(3) / 3]),
                  Edge(edges[2], lambda x: [(1 - x) / 2,
                                            np.sqrt(3) * (3 * x + 1) / 6], o=r),
                  Edge(edges[5], lambda x: [(- x - 1) / 2,
                                            np.sqrt(3) * (3 * -x + 1) / 6])])

face2 = Point(2, [Edge(edges[0], lambda x: [x, -np.sqrt(3) / 3]),
                  Edge(edges[4], lambda x: [(1 - x) / 2,
                                           np.sqrt(3) * (3 * x + 1) / 6]),
                  Edge(edges[3], lambda x: [(- x - 1) / 2,
                                           np.sqrt(3) * (3 * -x + 1) / 6])])
face3 = Point(2, [Edge(edges[0], lambda x: [x, -np.sqrt(3) / 3]),
                  Edge(edges[1], lambda x: [(1 - x) / 2,
                                           np.sqrt(3) * (3 * x + 1) / 6]),
                  Edge(edges[2], lambda x: [(- x - 1) / 2,
                                           np.sqrt(3) * (3 * -x + 1) / 6])])

face4 = Point(2, [Edge(edges[4], lambda x: [x, -np.sqrt(3) / 3]),
                  Edge(edges[5], lambda x: [(1 - x) / 2,
                                            np.sqrt(3) * (3 * x + 1) / 6], o=r),
                  Edge(edges[1], lambda x: [(- x - 1) / 2,
                                            np.sqrt(3) * (3 * -x + 1) / 6], o=r)])


lst3 = [Edge(face1, lambda x, y: [-0.5*x + (np.sqrt(3)/6)*y - 1/3, 0.5*x + (np.sqrt(3)/2)*y, -(np.sqrt(2)/2)*x + (np.sqrt(6)/6)*y + (np.sqrt(2) / 6)]),
        Edge(face2, lambda x, y: [x, -(np.sqrt(3)/3)*y - 1/3, (np.sqrt(6)/3)*y - (np.sqrt(2)/6)]),
        Edge(face3, lambda x, y: [x, (np.sqrt(3)/3)*y + 1/3,  (np.sqrt(6)/3)*y - (np.sqrt(2)/6)]),
        Edge(face4, lambda x, y: [-0.5*x - (np.sqrt(3)/6)*y + 1/3, -0.5*x + (np.sqrt(3)/2)*y, (np.sqrt(2)/2)*x + (np.sqrt(6)/6)*y + (np.sqrt(2) / 6)])]

tetrahedron = Point(3, lst3)

# print(tetrahedron.vertices())
# print(face1.cell_attachment(3))
# print(face1.cell_attachment(0))
# print(face1.cell_attachment(2))

# print(face2.cell_attachment(0))
# print(face2.cell_attachment(1))
# print(face2.cell_attachment(3))

# print(face3.cell_attachment(0))
# print(face3.cell_attachment(1))
# print(face3.cell_attachment(2))

# print(face4.cell_attachment(1))
# print(face4.cell_attachment(3))
# print(face4.cell_attachment(2))

# print(tetrahedron.cell_attachment(13)(twod[0][1], twod[0][2]))
# print(tetrahedron.cell_attachment(13)(twod[1][1], twod[1][2]))
# print(tetrahedron.cell_attachment(13)(twod[2][1], twod[2][2]))


# print(tetrahedron.cell_attachment(0))
# print(tetrahedron.cell_attachment(1))
# print(tetrahedron.cell_attachment(2))
# print(tetrahedron.cell_attachment(3))
# print(tetrahedron.cell_attachment(4)(-1))
# print(tetrahedron.cell_attachment(4)(1))
# print(tetrahedron.cell_attachment(5)(-1))
# print(tetrahedron.cell_attachment(5)(1))
# print(tetrahedron.cell_attachment(6)(-1))
# print(tetrahedron.cell_attachment(6)(1))
# print(tetrahedron.cell_attachment(7)(-1))
# print(tetrahedron.cell_attachment(7)(1))
# print(tetrahedron.cell_attachment(8)(-1))
# print(tetrahedron.cell_attachment(8)(1))
# print(tetrahedron.dimension)
# tetrahedron.hasse_diagram()
# A4 = GroupRepresentation(AlternatingGroup(4))
# s2tet = (S4/S2).add_cell(tetrahedron)
# for m in s2tet.members():
#     print(tetrahedron.permute_entities(m, 1))
#     print(tetrahedron.permute_entities(m, 2))

# breakpoint()
# xs = [lambda g: DOF(DeltaPairing(), PointKernel(g((0.3, 0.3, 0.3))))]
# dg1 = ElementTriple(tetrahedron, (P1, CellL2, "C0"),
#                     DOFGenerator(xs, S4, S1))
# ls = dg1.generate()
# print("num dofs ", dg1.num_dofs())
# fig = plt.figure()
# ax = fig.add_subplot(projection='3d')
# tetrahedron.plot3d(show=False, ax=ax)
# for dof in ls:
#     print(dof)
#     plotted = dof(lambda x, y, z: (x, y, z))
#     ax.scatter(plotted[0], plotted[1], plotted[2])
# plt.show()

# xs = [lambda g: DOF(DeltaPairing(), PointKernel(g((-1, 0, -1/np.sqrt(2)))))]
# dg1 = ElementTriple(tetrahedron, (P1, CellL2, "C0"),
#                     DOFGenerator(xs, S4/S2, S1))
# ls = dg1.generate()
# print("num dofs ", dg1.num_dofs())
# fig = plt.figure()
# ax = fig.add_subplot(projection='3d')
# tetrahedron.plot3d(show=False, ax=ax)
# for dof in ls:
#     print(dof)
#     plotted = dof(lambda x, y, z: (x, y, z))
#     ax.scatter(plotted[0], plotted[1], plotted[2], color="b")
# plt.show()

Z4 = GroupRepresentation(CyclicGroup(4))
print(Z4.base_group.generators)

Z3 = GroupRepresentation(CyclicGroup(3))
print(Z3.base_group.generators)

Z2 = GroupRepresentation(CyclicGroup(2))
print(Z2.base_group.generators)

D2 = GroupRepresentation(DihedralGroup(2))
print(D2.base_group.generators)

# # C6 = GroupRepresentation(CyclicGroup(6))

# A4 = GroupRepresentation(AlternatingGroup(4))
# print(A4.base_group.generators)

# a = Permutation(0, 1)(2, 3)
# b = Permutation(0, 2)(1, 3)
# c = Permutation(0, 3)(1, 2)
# K4 = GroupRepresentation(PermutationGroup(a, b, c))

# # breakpoint()
# print((D2 * C3).base_group.generators)

identity = MyTestFunction(lambda *x: x)
# # xs = [lambda g: DOF(DeltaPairing(), PointKernel(g((-0.6, 0, -1/np.sqrt(2)))))]

# # dg1 = ElementTriple(tetrahedron, (P1, CellL2, "C0"),
# #                     DOFGenerator(xs, A4, S1))
# # ls = dg1.generate()
# # print("num dofs ", dg1.num_dofs())
# # fig = plt.figure()
# # ax = fig.add_subplot(projection='3d')
# # tetrahedron.plot3d(show=False, ax=ax)
# # for dof in ls:
# #     print(dof)
# #     plotted = dof(identity)
# #     ax.scatter(plotted[0], plotted[1], plotted[2], color="b")
# # plt.show()

print("CG3 on tetrahedron")

xs = [lambda g: DOF(DeltaPairing(), PointKernel(g(())))]
dg0 = ElementTriple(vertices[0], (P0, CellL2, "C0"),
                    DOFGenerator(xs, S1, S1))

xs = [lambda g: DOF(DeltaPairing(), PointKernel(g((-0.3,)))),
      lambda g: DOF(DeltaPairing(), PointKernel(g((0.5,))))]
dg1_int = ElementTriple(edges[0], (P0, CellL2, "C0"),
                        DOFGenerator(xs, S1, S1))

xs = [lambda g: DOF(DeltaPairing(), PointKernel(g((0, 0))))]
dg0_face = ElementTriple(face1, (P0, CellL2, "C0"),
                        DOFGenerator(xs, S1, S1))


v_xs = [lambda g: immerse(g, tetrahedron, dg0, CellH1)]
verts = DOFGenerator(v_xs, S4 / S2, S1)

e_xs = [lambda g: immerse(g, tetrahedron, dg1_int, CellH1)]
edges = DOFGenerator(e_xs, Z3 * D2, S1)

f_xs = [lambda g: immerse(g, tetrahedron, dg0_face, CellH1)]
faces = DOFGenerator(f_xs, S4 / S2, S1)


cg3 = ElementTriple(tetrahedron, (P1, CellH1, "C0"),
                    [verts, edges, faces])
ls = cg3.generate()
fig = plt.figure()
ax = fig.add_subplot(projection='3d')

tetrahedron.plot3d(show=False, ax=ax)
count = 0
for dof in ls:
    print(dof)
    count += 1
    plotted = dof(identity)
    if dof.trace_entity.dimension == 1:
        color = "r"
    elif dof.trace_entity.dimension == 2:
        color = "g"
    else:
        color = "b"
    
    ax.scatter(plotted[0], plotted[1], plotted[2], color=color)
print("num dofs ", cg3.num_dofs())
print("actual num dofs ", count)
plt.show()