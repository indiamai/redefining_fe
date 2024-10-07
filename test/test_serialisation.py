from redefining_fe import *
from redefining_fe.serialisation import FETripleEncoder
from test_convert_to_fiat import create_cg1
import json

vert = Point(0)
edge = Point(1, [Point(0), Point(0)], vertex_num=2)
tri = n_sided_polygon(3)


def test_repeated_objs():
    repeated_edge = Point(1, [vert, vert], vertex_num=2)
    encoded = json.dumps(repeated_edge, cls=FETripleEncoder, indent=2)
    print(encoded)
    FETripleEncoder().encode(repeated_edge)


def test_serialisation():
    cells = [vert, edge, tri]

    for cell in cells:
        triple = create_cg1(cell)
        encoded = json.dumps(triple, cls=FETripleEncoder, indent=2)
        print(type(encoded))


def test_brackets():
    json_obj = FETripleEncoder()

    bracket_str = "start{bracket{ content { filler } more }}end"

    s, e = json_obj.bracket_matching(bracket_str)

    print(bracket_str[s:e])
