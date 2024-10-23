from redefining_fe import *
from redefining_fe.serialisation import FETripleEncoder, FETripleDecoder, bracket_matching, ElementSerialiser
from test_convert_to_fiat import create_cg1
import json

vert = Point(0)
edge = Point(1, [Point(0), Point(0)], vertex_num=2)
tri = n_sided_polygon(3)


# def test_simple():
#     print(vert)
#     encoded = json.dumps(vert, cls=FETripleEncoder, check_circular=False)
#     decoded = json.loads(encoded, cls=FETripleDecoder)
#     xs = [DOF(DeltaPairing(), PointKernel(()))]
#     dg0 = ElementTriple(decoded, (P0, CellL2, C0), DOFGenerator(xs, S1, S1))
#     for dof in dg0.generate():
#         print(dof)
#     print(decoded)

#     encoded = json.dumps(dg0, cls=FETripleEncoder)
    # decoded = json.loads(encoded, cls=FETripleDecoder)

    # encoded = json.dumps(edge, cls=FETripleEncoder, check_circular=False)
    # decoded = json.loads(encoded, cls=FETripleDecoder)
    # print(decoded)
    # xs = [DOF(DeltaPairing(), PointKernel((-1,)))]
    # dg1 = ElementTriple(decoded, (P1, CellL2, C0), DOFGenerator(xs, S2, S1))
    # for dof in dg1.generate():
    #     print(dof)


def test_repeated_objs():
    repeated_edge = Point(1, [vert, vert], vertex_num=2)
    repeated_tri = Point(2, [repeated_edge, repeated_edge, repeated_edge], vertex_num=3)
    converter = ElementSerialiser()
    encoded = converter.encode(repeated_edge)
    decoded = converter.decode(encoded)
    print(encoded)
    print(decoded)


def test_serialisation():
    cells = [vert, edge, tri]

    for cell in cells:
        triple = create_cg1(cell)
        converter = ElementSerialiser()
        encoded = converter.encode(triple)
        print(encoded)
        decoded = converter.decode(encoded)
        print(decoded)


def test_brackets():
    bracket_str = "start{bracket{ content { filler } more }}end"

    s, e = bracket_matching(bracket_str)

    assert (bracket_str[s:e] == "{bracket{ content { filler } more }}")
