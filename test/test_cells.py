from redefining_fe import *
import pytest
import numpy as np


@pytest.fixture(scope='module', params=[0, 1, 2])
def C(request):
    dim = request.param
    if dim == 0:
        return Point(0)
    elif dim == 1:
        return Point(1, [Point(0), Point(0)], vertex_num=2)
    elif dim == 2:
        return n_sided_polygon(3)


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


def test_orientation():
    cell = Point(1, [Point(0), Point(0)], vertex_num=2)
    print(cell.get_topology())
    for g in cell.group.members():
        if not g.perm.is_Identity:
            oriented = cell.orient(g)
            assert np.allclose(np.array(oriented.basis_vectors(return_coords=True)[0]), -1)


def test_sub_basis_vectors():
    cell = n_sided_polygon(3)

    edges = cell.edges(get_class=True)
    print(cell.vertices())
    print(cell.vertices(return_coords=True))
    for e in edges:
        print(e)
        print(e.vertices())
        print(cell.basis_vectors(entity=e))
