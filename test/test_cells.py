import numpy as np
from src.cell_complex.cells import Point, Edge
import pytest


@pytest.fixture(scope='module', params=[0, 1, 2])
def C(request):
    dim = request.param
    if dim == 0:
        return Point(0)
    elif dim == 1:
        vertices = [Point(0), Point(0)]
        return Point(1, [Edge(vertices[0], lambda: (-1,)),
                         Edge(vertices[1], lambda: (1,))])
    elif dim == 2:
        vertices = []
        for i in range(3):
            vertices.append(Point(0))
        edges = []
        edges.append(
            Point(1, [Edge(vertices[0], lambda: (-1,)),
                      Edge(vertices[1], lambda: (1,))]))
        edges.append(
            Point(1, [Edge(vertices[0], lambda: (1,)),
                      Edge(vertices[2], lambda: (-1,))]))
        edges.append(
            Point(1, [Edge(vertices[1], lambda: (-1,)),
                      Edge(vertices[2], lambda: (1,))]))

        return Point(2, [Edge(edges[0], lambda x: [x, -np.sqrt(3)/3]),
                         Edge(edges[1], lambda x: [(-x - 1)/2,
                                                   np.sqrt(3)*(3*-x + 1)/6]),
                         Edge(edges[2], lambda x: [(1 - x)/2,
                                                   np.sqrt(3)*(3*x + 1)/6])])


def test_vertices(C):
    verts = C.vertices()
    assert len(verts) == C.dimension + 1


def test_basis_vectors(C):
    if C.dimension == 0:
        with pytest.raises(ValueError):
            bv_ids = C.basis_vectors()
        with pytest.raises(ValueError):
            bv_coords = C.basis_vectors(return_coords=True)
    else:
        bv_ids = C.basis_vectors()
        bv_coords = C.basis_vectors(return_coords=True)
        assert len(bv_ids) == len(bv_coords)
