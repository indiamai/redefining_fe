from redefining_fe import *
import numpy as np
import matplotlib.pyplot as plt


triangle_nums = [sum(range(i)) for i in range(2, 20)]


def triangle_coords(n, verts=[(-1, -np.sqrt(3)/3), (0, 2*np.sqrt(3)/3), (1, -np.sqrt(3)/3)], scale=0.8):
    assert n in triangle_nums
    if n == 1:
        return [(0, 0)]
    on_edge = triangle_nums.index(n) - 1
    sub_verts = [(scale*x, scale*y) for (x, y) in verts]
    if on_edge > 0:
        ratio = 1 / (on_edge + 1)
        for i in range(on_edge):
            sub_verts += [(sub_verts[0][0] + ratio*(i+1)*(sub_verts[1][0] - sub_verts[0][0]), sub_verts[0][1] + ratio*(i+1)*(sub_verts[1][1] - sub_verts[0][1]))]
            sub_verts += [(sub_verts[0][0] + ratio*(i+1)*(sub_verts[2][0] - sub_verts[0][0]), sub_verts[0][1] + ratio*(i+1)*(sub_verts[2][1] - sub_verts[0][1]))]
            sub_verts += [(sub_verts[2][0] + ratio*(i+1)*(sub_verts[1][0] - sub_verts[2][0]), sub_verts[2][1] + ratio*(i+1)*(sub_verts[1][1] - sub_verts[2][1]))]
        if (n - len(sub_verts)) > 0:
            int_scale = (sub_verts[0][0] + ratio*3*(sub_verts[1][0] - sub_verts[0][0])) / sub_verts[0][0]
            interior = triangle_coords((n - len(sub_verts)), sub_verts[0:3], int_scale)
            sub_verts += interior
    assert len(sub_verts) == n
    return sub_verts


def CR_n(deg):
    cell = polygon(3)
    points = np.polynomial.legendre.leggauss(deg)[0]
    print(points)
    Pk = PolynomialSpace(3)
    sym_points = [DOF(DeltaPairing(), PointKernel((pt,))) for pt in points[:len(points)//2]]
    print(sym_points)
    if len(points) % 2 == 0:
        edge_dg0 = ElementTriple(cell.edges(get_class=True)[0], (Pk, CellL2, C0), [DOFGenerator(sym_points, S2, S1),
                                                                                   DOFGenerator([DOF(DeltaPairing(), PointKernel((0,)))], S1, S1)])
    else:
        edge_dg0 = ElementTriple(cell.edges(get_class=True)[0], (Pk, CellL2, C0), [DOFGenerator(sym_points, S2, S1)])
    edge_xs = [immerse(cell, edge_dg0, TrH1)]
    center = [DOF(DeltaPairing(), PointKernel((0, 0)))]

    return ElementTriple(cell, (Pk, CellL2, C0), [DOFGenerator(edge_xs, C3, S1), DOFGenerator(center, S1, S1)])


coords = triangle_coords(36)
edge_coords = [(-1, -np.sqrt(3)/3), (0, 2*np.sqrt(3)/3), (1, -np.sqrt(3)/3)]
print(triangle_nums)
ax = plt.gca()
ax.scatter([e[0] for e in coords], [e[1] for e in coords], color="black")
ax.scatter([e[0] for e in edge_coords], [e[1] for e in edge_coords], color="blue")
ax.figure.savefig("triangle.png")
CR_n(3)
