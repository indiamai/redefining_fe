from redefining_fe import *
from test_convert_to_fiat import create_cg1
import pytest

vert = Point(0)
edge = Point(1, [Point(0), Point(0)], vertex_num=2)
tri = n_sided_polygon(3)



def test_serialisation():
    cells = [vert, edge, tri]

    for cell in cells:
        triple = create_cg1(cell)
        triple.to_json(str(cell) + ".json")

    