from redefining_fe import *
from redefining_fe.serialisation import FETripleEncoder, bracket_matching
# FETripleDecoder,
from test_convert_to_fiat import create_cg1
import json

vert = Point(0)
edge = Point(1, [Point(0), Point(0)], vertex_num=2)
tri = n_sided_polygon(3)


def test_repeated_objs():
    repeated_edge = Point(1, [vert, vert], vertex_num=2)
    repeated_tri = Point(2, [repeated_edge, repeated_edge, repeated_edge], vertex_num=3)
    encoded = json.dumps(repeated_tri, cls=FETripleEncoder, indent=2, check_circular=False)
    print(encoded)


def test_serialisation():
    cells = [vert, edge, tri]

    for cell in cells:
        triple = create_cg1(cell)
        encoded = json.dumps(triple, cls=FETripleEncoder, indent=2, check_circular=False)
        print(encoded)
        # decoded = json.loads(encoded, cls=FETripleDecoder)
        # print(type(decoded))


def test_brackets():
    bracket_str = "start{bracket{ content { filler } more }}end"

    s, e = bracket_matching(bracket_str)

    assert (bracket_str[s:e] == "{bracket{ content { filler } more }}")
