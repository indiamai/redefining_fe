# Examples of Elements in 3 dimensions
from firedrake import *
import numpy as np
from redefining_fe import *

vertices = []
for i in range(4):
    vertices.append(Point(0))
edges = []
edges.append(
    Point(1, vertex_num=2, edges=[vertices[0], vertices[1]]))
edges.append(
    Point(1, vertex_num=2, edges=[vertices[1], vertices[2]]))
edges.append(
    Point(1, vertex_num=2, edges=[vertices[2], vertices[0]]))
edges.append(
    Point(1, vertex_num=2, edges=[vertices[3], vertices[0]]))
edges.append(
    Point(1, vertex_num=2, edges=[vertices[1], vertices[3]]))
edges.append(
    Point(1, vertex_num=2, edges=[vertices[2], vertices[3]]))


face1 = Point(2, vertex_num=3, edges=[edges[5], edges[3], edges[2]], edge_orientations={2: r})
face2 = Point(2, vertex_num=3, edges=[edges[3], edges[0], edges[4]])
face3 = Point(2, vertex_num=3, edges=[edges[2], edges[0], edges[1]])
face4 = Point(2, vertex_num=3, edges=[edges[1], edges[4], edges[5]], edge_orientations={0: r, 2: r})

tetrahedron = Point(3, vertex_num=4, edges=[face3, face1, face4, face2])
tetra = tetrahedron
# tetrahedron.plot3d()
print(tetrahedron.vertices())

xs = [DOF(DeltaPairing(), PointKernel((0.3, 0.3, 0.3)))]
dg0 = ElementTriple(tetra, (P1, CellL2, "C0"),
                    DOFGenerator(xs, Z4, S1))
ls = dg0.generate()
# dg0.plot()
# print("num dofs ", dg0.num_dofs())
# for dof in ls:
#     print(dof)


xs = [DOF(DeltaPairing(), PointKernel(tuple(tetra.vertices(return_coords=True)[0])))]
dg1 = ElementTriple(tetra, (P1, CellL2, "C0"),
                    DOFGenerator(xs, Z4, S1))
ls = dg1.generate()
# dg1.plot()
# print("num dofs ", dg1.num_dofs())
# for dof in ls:
#     print(dof)


xs = [DOF(DeltaPairing(), PointKernel((-1/np.sqrt(2) + 0.4, -1/np.sqrt(2)+0.4, 1/np.sqrt(2) - 0.2)))]

dg1 = ElementTriple(tetra, (P1, CellL2, "C0"),
                    DOFGenerator(xs, A4, S1))
ls = dg1.generate()
# dg1.plot()
# for dof in ls:
#     print(dof)
# print("num dofs ", dg1.num_dofs())


print("CG3 on tetrahedron")

xs = [DOF(DeltaPairing(), PointKernel(()))]
dg0 = ElementTriple(vertices[0], (P0, CellL2, "C0"),
                    DOFGenerator(xs, S1, S1))

xs = [DOF(DeltaPairing(), PointKernel((-1/3,)))]
dg1_int = ElementTriple(edges[0], (P0, CellL2, "C0"),
                        DOFGenerator(xs, S2, S1))

xs = [DOF(DeltaPairing(), PointKernel((0, 0)))]
dg0_face = ElementTriple(face1, (P0, CellL2, "C0"),
                         DOFGenerator(xs, S1, S1))

v_xs = [immerse(tetra, dg0, CellH1)]
cgverts = DOFGenerator(v_xs, Z4, S1)

e_xs = [immerse(tetra, dg1_int, CellH1)]
cgedges = DOFGenerator(e_xs, A4, S1)

f_xs = [immerse(tetra, dg0_face, CellH1)]
cgfaces = DOFGenerator(f_xs, S4, S1)


cg3 = ElementTriple(tetra, (P1, CellH1, "C0"),
                    [cgverts, cgedges, cgfaces])

ls = cg3.generate()
# cg3.plot()
# for dof in ls:
#     print(dof)
# print("num dofs ", cg3.num_dofs())


print("Raviart Thomas")
xs = [DOF(L2InnerProd(), PolynomialKernel(lambda x: 1))]
dofs = DOFGenerator(xs, S1, S2)
face_vec = ElementTriple(face1, (P1, CellHDiv, "C0"), dofs)
ls = face_vec.generate()

im_xs = [immerse(tetra, face_vec, CellHDiv)]
face = DOFGenerator(im_xs, S4, S4)

rt1 = ElementTriple(tetra, (P1, CellHDiv, "C0"),
                    [face])
ls = rt1.generate()
# rt1.plot()
# for dof in ls:
#     print(dof)
# print("num dofs ", rt1.num_dofs())

print("Nedelec")
xs = [DOF(L2InnerProd(), PolynomialKernel(lambda x: 1))]
dofs = DOFGenerator(xs, S1, S2)
int_ned = ElementTriple(edges[0], (P1, CellHCurl, "C0"), dofs)
ls = int_ned.generate()

im_xs = [immerse(tetra, int_ned, CellHCurl)]
edge = DOFGenerator(im_xs, A4, Z4)

ned = ElementTriple(tetra, (P1, CellHCurl, "C0"),
                    [edge])
ls = ned.generate()
ned.plot()
for dof in ls:
    print(dof)
print("num dofs ", ned.num_dofs())
