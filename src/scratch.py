import numpy as np
from cell_complex.cells import Point, Edge
import matplotlib.pyplot as plt
import sympy as sp
from groups.new_groups import r, rot, S1, S2, S3, D4, C4
from dof_lang.dof import DeltaPairing, DOF, L2InnerProd, MyTestFunction, PointKernel
from triples import ElementTriple, DOFGenerator, immerse
from spaces.element_sobolev_spaces import CellH1, CellL2, CellHDiv, CellHCurl, CellH2, CellH3
from spaces.polynomial_spaces import P0, P1, P2, P3, Q2, VectorPolynomialSpace

def compute_scaled_verts(d, n):
    if d == 2:
        source = np.array([0, 1])
        rot_coords = [source for i in range(0, n)]

        rot_mat = np.array([[np.cos((2*np.pi)/n), -np.sin((2*np.pi)/n)],[np.sin((2*np.pi)/n), np.cos((2*np.pi)/n)]])
        for i in range(1, n):
            rot_coords[i] = np.matmul(rot_mat, rot_coords[i-1])
        xdiff, ydiff = (rot_coords[0][0] - rot_coords[1][0],
                        rot_coords[0][1] - rot_coords[1][1])
        scale = 2 / np.sqrt(xdiff**2 + ydiff**2)
        scaled_coords = np.array([[scale*x, scale*y] for (x, y) in rot_coords])
        return scaled_coords
    elif d == 3:
        if n == 4:
            A = [-1, 1, -1]
            B = [1, -1, -1]
            C = [1, 1, 1]
            D = [-1, -1, 1]
            coords = [A, B, C, D]
            face1 = np.array([A, D, C])
            face2 = np.array([A, B, D])
            face3 = np.array([A, C, B])
            face4 = np.array([B, D, C])
            faces = [face1, face2, face3, face4]
        elif n == 8:
            coords = []
            faces = [[] for i in range(6)]
            for i in [-1, 1]:
                for j in [-1, 1]:
                    for k in [-1, 1]:
                        coords.append([i, j, k])
            
            for j in [-1, 1]:
                for k in [-1, 1]:
                    faces[0].append([1, j, k])
                    faces[1].append([-1, j, k])
                    faces[2].append([j, 1, k])
                    faces[3].append([j, -1, k])
                    faces[4].append([j, k, 1])
                    faces[5].append([j, k, -1])

        else:
            raise ValueError("Polyhedron with {} vertices not supported".format(n))

        # xdiff, ydiff, zdiff = (coords[0][0] - coords[1][0],
        #                        coords[0][1] - coords[1][1],
        #                        coords[0][2] - coords[1][2])
        # scale = 2 / np.sqrt(xdiff**2 + ydiff**2 + zdiff**2)
        # scaled_coords = np.array([[scale*x, scale*y, scale*z] for (x, y, z) in coords])
        # scaled_faces = np.array([[[scale*x, scale*y, scale*z] for (x, y, z) in face] for face in faces])

        # return scaled_coords, scaled_faces
        return coords, faces
    else:
        raise ValueError("Dimension {} not supported".format(d))

def n_sided_polygon(n):
    vertices = []
    for i in range(n):
        vertices.append(Point(0))
    edges = []
    for i in range(n):
        edges.append(
            Point(1, compute_attachments(2, 1, [vertices[(i+1) % n],
                                                vertices[(i+2) % n]])))

    return Point(2, compute_attachments(n, 2, edges))

def compute_attachments(d, n, points, orientations={}):
    if d == 1:
        edges = [Edge(points[0], lambda: (-1,)),
                    Edge(points[1], lambda: (1,))]
    if d == 2:
        coords = compute_scaled_verts(2, n)
        edges = []

        for i in range(n):
            a, b = coords[i]
            c, d = coords[(i + 1) % n]

            if i in orientations.keys():
                edges.append(Edge(points[i], construct_attach_2d(a, b, c, d), o=orientations[i]))
            else:
                edges.append(Edge(points[i], construct_attach_2d(a, b, c, d)))
    if d == 3:
        coords, faces = compute_scaled_verts(3, n)
        
        coords_2d = np.c_[np.ones(len(faces[0])), compute_scaled_verts(2, len(faces[0]))]
        
        res = []
        edges = []

        for i in range(len(faces)):
            res = np.linalg.solve(coords_2d, faces[i])

            res_fn = construct_attach_3d(res)
            assert np.allclose(res_fn(coords_2d[0][1], coords_2d[0][2]), faces[i][0])
            assert np.allclose(res_fn(coords_2d[1][1], coords_2d[1][2]), faces[i][1])
            assert np.allclose(res_fn(coords_2d[2][1], coords_2d[2][2]), faces[i][2])
            if i in orientations.keys():
                edges.append(Edge(points[i], construct_attach_3d(res), o=orientations[i]))
            else:
                edges.append(Edge(points[i], construct_attach_3d(res)))

        # breakpoint()
    return edges



def construct_attach_2d(a, b, c, d):
    return lambda x: [((c-a)/2)*(x+1) + a, ((d-b)/2)*(x+1) + b]


def construct_attach_3d(res):
    x = sp.Symbol("x")
    y = sp.Symbol("y")
    xy = sp.Matrix([1, x, y])
    # print(res.T*xy)
    # breakpoint()
    return lambda x, y: np.array((xy.T * res).subs({"x": x, "y": y}), dtype=float)




def r(x):
    # reflection in the first component
    if isinstance(x, int):
        x = [x]
    x_list = list(x)
    x_list[0] = -x_list[0]
    return tuple(x_list)

# tri = n_sided_polygon(3)
# # tri2.plot()
# print(tri.vertices(return_coords=True))
# print(tri.group.size())

# square = n_sided_polygon(4)
# # tri2.plot()
# print(square.vertices(return_coords=True))
# print(square.group.size())

# xs = [DOF(DeltaPairing(), PointKernel((-np.sqrt(2), 0)))]

# dg1 = ElementTriple(square, (P1, CellL2, "C0"),
#                     DOFGenerator(xs, C4, S1))
# ls = dg1.generate()
# for l in ls:
#     print(l)
# dg1.plot()

# v_xs = [immerse(square, dg0, CellH1)]
# v_dofs = DOFGenerator(v_xs, S3/S2, S1)

# e_xs = [immerse(a4, dg0_int, CellH1)]
# e_dofs = DOFGenerator(e_xs, S3, S1)

# i_xs = [lambda g: DOF(DeltaPairing(), PointKernel(g((0, 0))))]
# i_dofs = DOFGenerator(i_xs, S1, S1)

# cg3 = ElementTriple(a4, (P3, CellH1, "C0"),
#                     [v_dofs, e_dofs, i_dofs])

# phi_0 = MyTestFunction(lambda x, y: (x, y))
# ls = cg3.generate()
# print("num dofs ", cg3.num_dofs())
# for dof in ls:
#     print(dof)
#     print(dof.eval(phi_0))
# cg3.plot()
# hex = n_sided_polygon(6)
# hex.plot()
# print(hex.vertices(return_coords=True))

# print(compute_scaled_verts_and_faces_3d(4))
# print(compute_scaled_verts_and_faces_3d(8))
# print(np.sqrt(2)/2)
# lst = compute_attachments(4, 3, ["A","B", "C", "D"])
# for a in lst:
#     print(a(-1, -np.sqrt(2)/2))


vertices1 = []
for i in range(4):
    vertices1.append(Point(0))
edges1 = []
edges1.append(
    Point(1, vertex_num=2, edges=[vertices1[0], vertices1[1]]))  # 4
edges1.append(
    Point(1, vertex_num=2, edges=[vertices1[1], vertices1[2]]))  # 5
edges1.append(
    Point(1, vertex_num=2, edges=[vertices1[2], vertices1[0]]))  # 6
edges1.append(
    Point(1, vertex_num=2, edges=[vertices1[3], vertices1[0]]))  # 7
edges1.append(
    Point(1, vertex_num=2, edges=[vertices1[1], vertices1[3]]))  # 8
edges1.append(
    Point(1, vertex_num=2, edges=[vertices1[2], vertices1[3]]))  # 9
# aa = Point(2, vertex_num=3, edges=[edges1[1], edges1[4], edges1[5]], edge_orientations= {0:r, 2: r})
# aa = Point(2, vertex_num=3, edges=[edges1[2], edges1[0], edges1[1]])
# aa =  Point(2, vertex_num=3, edges=[edges1[3], edges1[0], edges1[4]])

# aa = Point(2, vertex_num=3, edges=[edges1[5], edges1[3], edges1[2]], edge_orientations={0: r, 1:r})

# for a in aa.d_entities(3):
#     print(a)
#     print(aa.cell_attachment(a)())

face10 = Point(2, vertex_num=3, edges=[edges1[5], edges1[3], edges1[2]], edge_orientations={2: r})
face11 = Point(2, vertex_num=3, edges=[edges1[3], edges1[0], edges1[4]])
face12 = Point(2, vertex_num=3, edges=[edges1[2], edges1[0], edges1[1]])
face13 = Point(2, vertex_num=3, edges=[edges1[1], edges1[4], edges1[5]], edge_orientations= {0:r, 2: r})
# e10 = compute_attachments(3, 2, [edges1[5], edges1[3], edges1[2]])
# e11 = compute_attachments(3, 2, [edges1[3], edges1[0], edges1[4]])
# e12 = compute_attachments(3, 2, [edges1[2], edges1[0], edges1[1]])
# e13 = compute_attachments(3, 2, [edges1[1], edges1[4], edges1[5]])

# # breakpoint()
# face10 = Point(2, compute_attachments(3, 2, [edges1[5], edges1[3], edges1[2]]))
# face11 = Point(2, compute_attachments(3, 2, [edges1[3], edges1[0], edges1[4]]))
# face12 = Point(2, compute_attachments(3, 2, [edges1[2], edges1[0], edges1[1]]))
# face13 = Point(2, compute_attachments(3, 2, [edges1[1], edges1[4], edges1[5]]))

f12, f10, f13, f11 = compute_attachments(3, 4, [face12, face10, face13, face11])
tetra = Point(3, vertex_num=4, edges=[face12, face10, face13, face11])
# # tetra.hasse_diagram()
# verts = compute_scaled_verts(2, 3)
# print(verts)

# print("vert 0")
# print("12", f12(*verts[1]))
# print("11", f11(*verts[1]))
# print("10", f10(*verts[2]))

# print("vert 1")
# print("13", f13(*verts[1]))
# print("12", f12(*verts[2]))
# print("11", f11(*verts[2]))

# print("vert 2")
# print("13", f13(*verts[0]))
# print("12", f12(*verts[0]))
# print("10", f10(*verts[0]))

# print("vert 3")
# print("13", f13(*verts[2]))
# print("11", f11(*verts[0]))
# print("10", f10(*verts[1]))
# breakpoint()

# print(tetra.basis_vectors(return_coords=True))
# face1 = Point(2, [Edge(edges[3], lambda x: [x, -np.sqrt(3) / 3]),
#                   Edge(edges[2], lambda x: [(1 - x) / 2,
#                                             np.sqrt(3) * (3 * x + 1) / 6], o=r),
#                   Edge(edges[5], lambda x: [(- x - 1) / 2,
#                                             np.sqrt(3) * (3 * -x + 1) / 6])])

# face2 = Point(2, [Edge(edges[0], lambda x: [x, -np.sqrt(3) / 3]),
#                   Edge(edges[4], lambda x: [(1 - x) / 2,
#                                            np.sqrt(3) * (3 * x + 1) / 6]),
#                   Edge(edges[3], lambda x: [(- x - 1) / 2,
#                                            np.sqrt(3) * (3 * -x + 1) / 6])])
# face3 = Point(2, [Edge(edges[0], lambda x: [x, -np.sqrt(3) / 3]),
#                   Edge(edges[1], lambda x: [(1 - x) / 2,
#                                            np.sqrt(3) * (3 * x + 1) / 6]),
#                   Edge(edges[2], lambda x: [(- x - 1) / 2,
#                                            np.sqrt(3) * (3 * -x + 1) / 6])])

# face4 = Point(2, [Edge(edges[4], lambda x: [x, -np.sqrt(3) / 3]),
#                   Edge(edges[5], lambda x: [(1 - x) / 2,
#                                             np.sqrt(3) * (3 * x + 1) / 6], o=r),
#                   Edge(edges[1], lambda x: [(- x - 1) / 2,
#                                             np.sqrt(3) * (3 * -x + 1) / 6], o=r)])
