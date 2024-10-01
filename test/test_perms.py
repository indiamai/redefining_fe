from redefining_fe import *
from test_convert_to_fiat import create_cg1, create_dg1, create_cg2
from test_2d_examples_docs import construct_cg3
import pytest

vert = Point(0)
edge = Point(1, [Point(0), Point(0)], vertex_num=2)
tri = n_sided_polygon(3)


@pytest.mark.parametrize("cell", [edge])
def test_cg_perms(cell):
    cg1 = create_cg1(cell)
    cg1.make_basix_style_perms()


@pytest.mark.parametrize("cell", [edge])
def test_cg2_perms(cell):
    cg2 = create_cg2(cell)
    cg2.make_basix_style_perms()


@pytest.mark.parametrize("cell", [tri])
def test_cg3_perms(cell):
    cg3 = construct_cg3(cell)
    cg3.make_basix_style_perms()


@pytest.mark.parametrize("cell", [vert, edge])
def test_dg_perms(cell):
    cg1 = create_dg1(cell)
    cg1.make_basix_style_perms()
