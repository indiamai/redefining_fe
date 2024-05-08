# Examples of Elements in 3 dimensions
from firedrake import *
import numpy as np
from sympy.combinatorics import PermutationGroup, Permutation
from sympy.combinatorics.named_groups import SymmetricGroup, DihedralGroup, CyclicGroup, AlternatingGroup
from groups.new_groups import r, r_y, rot, S1, S2, S3, D4, C4, S4
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
twod = np.array([[1, -1, -np.sqrt(3)/3],
                 [1, 1, -np.sqrt(3)/3],
                 [1, 0, 2*np.sqrt(3)/3],])

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
# tetrahedron.plot3d()
# breakpoint()

s4 = SymmetricGroup(4)
print(s4.elements)
print(s4.generators)


print(S4.generators)
print(S4.members)
breakpoint()