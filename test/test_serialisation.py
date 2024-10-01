from redefining_fe import *
from redefining_fe.serialisation import FETripleEncoder
from test_convert_to_fiat import create_cg1
import pytest
import json

vert = Point(0)
edge = Point(1, [Point(0), Point(0)], vertex_num=2)
tri = n_sided_polygon(3)

def test_repeated_objs():
    repeated_edge = Point(1, [vert, vert], vertex_num=2)
    encoded = json.dumps(repeated_edge, cls=FETripleEncoder)
    print(encoded)

def test_serialisation():
    cells = [vert, edge, tri]

    for cell in cells:
        triple = create_cg1(cell)
        print(triple.__dict__)
        # triple.to_json(str(cell) + ".json")
        encoded = json.dumps(vert, cls=FETripleEncoder)
        print(encoded)