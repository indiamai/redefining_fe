import numpy as np
from cell_complex.cells import Point, Edge
import matplotlib.pyplot as plt
import sympy as sp

def compute_scaled_verts_2d(n):
    source = np.array([0, 1])
    rot_coords = [source for i in range(0, n)]

    rot_mat = np.array([[np.cos((2*np.pi)/n), -np.sin((2*np.pi)/n)],[np.sin((2*np.pi)/n), np.cos((2*np.pi)/n)]])
    for i in range(1, n):
        rot_coords[i] = np.matmul(rot_mat, rot_coords[i-1])
    xdiff, ydiff = (rot_coords[0][0] - rot_coords[1][0], rot_coords[0][1] - rot_coords[1][1])
    scale = 2 / np.sqrt(xdiff**2 + ydiff**2)
    scaled_coords = np.array([[scale*x, scale*y] for (x, y) in rot_coords])
    return scaled_coords


def compute_scaled_verts_and_faces_3d(n):
    if n == 4:
        A = [-1, 1, -1]
        B = [1, -1, -1]
        C = [1, 1, 1]
        D = [-1, -1, 1]
        coords = [A, B, C, D]
        face1 = np.array([D, A, C])
        face2 = np.array([A, B, D])
        face3 = np.array([A, B, C])
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

    xdiff, ydiff, zdiff = (coords[0][0] - coords[1][0],
                           coords[0][1] - coords[1][1],
                           coords[0][2] - coords[1][2])
    scale = 2 / np.sqrt(xdiff**2 + ydiff**2 + zdiff**2)
    scaled_coords = np.array([[scale*x, scale*y, scale*z] for (x, y, z) in coords])
    scaled_faces = np.array([[[scale*x, scale*y, scale*z] for (x, y, z) in face] for face in faces])

    return scaled_coords, scaled_faces

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

def compute_attachments(n, dimension, points, orientations={}):
    if dimension == 1:
        edges = [Edge(points[0], lambda: (-1,)),
                    Edge(points[1], lambda: (1,))]
    if dimension == 2:
        coords = compute_scaled_verts_2d(n)
        edges = []

        for i in range(n):
            a, b = coords[i]
            c, d = coords[(i + 1) % n]

            if i in orientations.keys():
                edges.append(Edge(points[i], construct_attach_2d(a, b, c, d), o=orientations[i]))
            else:
                edges.append(Edge(points[i], construct_attach_2d(a, b, c, d)))
    if dimension == 3:
        
        coords, faces = compute_scaled_verts_and_faces_3d(n)
        coords_2d = np.c_[np.ones(len(faces[0])), compute_scaled_verts_2d(len(faces[0]))]
        
        res = []
        edges = []

        for i in range(len(faces)):
            res = np.linalg.solve(coords_2d.T, faces[i])
            if i in orientations.keys():
                edges.append(Edge(points[i], construct_attach_3d(res), o=orientations[i]))
            else:
                edges.append(Edge(points[i], construct_attach_3d(res)))
    return edges


def construct_attach_2d(a, b, c, d):
    return lambda x: [((c-a)/2)*(x+1) + a, ((d-b)/2)*(x+1) + b]

def construct_attach_3d(res):
    x = sp.Symbol("x")
    y = sp.Symbol("y")
    xy = sp.Matrix([1, x, y])
    print(res.T*xy)
    return lambda x, y: list(np.array((res*xy).subs({"x": x, "y": y}), dtype=float))


def r(x):
    # reflection in the first component
    if isinstance(x, int):
        x = [x]
    x_list = list(x)
    x_list[0] = -x_list[0]
    return tuple(x_list)


# vertices1 = []
# for i in range(3):
#     vertices1.append(Point(0))
# edges1 = []
# edges1.append(
#     Point(1, compute_attachments(2, 1, [vertices1[1], vertices1[2]])))
# edges1.append(
#     Point(1, compute_attachments(2, 1, [vertices1[2], vertices1[0]])))
# edges1.append(
#     Point(1, compute_attachments(2, 1, [vertices1[0], vertices1[1]])))

# new_tri = Point(2, compute_attachments(3, 2, edges1))


# new_tri.plot()

# tri2 = n_sided_polygon(3)
# # tri2.plot()
# print(tri2.vertices(return_coords=True))

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
    Point(1, compute_attachments(2, 1, [vertices1[0], vertices1[1]])))  # 4
edges1.append(
    Point(1, compute_attachments(2, 1, [vertices1[1], vertices1[2]])))  # 5
edges1.append(
    Point(1, compute_attachments(2, 1, [vertices1[2], vertices1[0]])))  # 6
edges1.append(
    Point(1, compute_attachments(2, 1, [vertices1[3], vertices1[0]])))  # 7
edges1.append(
    Point(1, compute_attachments(2, 1, [vertices1[1], vertices1[3]])))  # 8
edges1.append(
    Point(1, compute_attachments(2, 1, [vertices1[2], vertices1[3]]))) # 9
aa = compute_attachments(3, 2, [edges1[1], edges1[4], edges1[5]])

for a in aa:
    print(a(-1))
    print(a(1))

face10 = Point(2, compute_attachments(3, 2, [edges1[5], edges1[3], edges1[2]], {2: r}))
face11 = Point(2, compute_attachments(3, 2, [edges1[3], edges1[0], edges1[4]]))
face12 = Point(2, compute_attachments(3, 2, [edges1[2], edges1[0], edges1[1]]))
face13 = Point(2, compute_attachments(3, 2, [edges1[1], edges1[4], edges1[5]], {0:r, 2: r}))
e10 = compute_attachments(3, 2, [edges1[5], edges1[3], edges1[2]])
e11 = compute_attachments(3, 2, [edges1[3], edges1[0], edges1[4]])
e12 = compute_attachments(3, 2, [edges1[2], edges1[0], edges1[1]])
e13 = compute_attachments(3, 2, [edges1[1], edges1[4], edges1[5]])

# # breakpoint()
# face10 = Point(2, compute_attachments(3, 2, [edges1[5], edges1[3], edges1[2]]))
# face11 = Point(2, compute_attachments(3, 2, [edges1[3], edges1[0], edges1[4]]))
# face12 = Point(2, compute_attachments(3, 2, [edges1[2], edges1[0], edges1[1]]))
# face13 = Point(2, compute_attachments(3, 2, [edges1[1], edges1[4], edges1[5]]))

f10, f11, f12, f13 = compute_attachments(4, 3, [face10, face11, face12, face12])
tetra = Point(3, compute_attachments(4, 3, [face13, face11, face12, face10])) 
# tetra.hasse_diagram()
verts = compute_scaled_verts_2d(3)

print("f10")
print(verts[0], f10(*verts[0]))
print(verts[1], f10(*verts[1]))
print(verts[2], f10(*verts[2]))

print("f11")
print(verts[0], f11(*verts[0]))
print(verts[1], f11(*verts[1]))
print(verts[2], f11(*verts[2]))

print("f12")
print(verts[0], f12(*verts[0]))
print(verts[1], f12(*verts[1]))
print(verts[2], f12(*verts[2]))

print("f13")
print(verts[0], f13(*verts[0]))
print(verts[1], f13(*verts[1]))
print(verts[2], f13(*verts[2]))
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
