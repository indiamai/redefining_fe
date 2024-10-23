from redefining_fe import *
from redefining_fe.serialisation import bracket_matching, ElementSerialiser
from test_convert_to_fiat import create_cg1
import numpy as np

vert = Point(0)
edge = Point(1, [Point(0), Point(0)], vertex_num=2)
tri = n_sided_polygon(3)


def test_simple():
    converter = ElementSerialiser()
    encoded = converter.encode(vert)
    decoded = converter.decode(encoded)
    xs = [DOF(DeltaPairing(), PointKernel(()))]
    dg0 = ElementTriple(decoded, (P0, CellL2, C0), DOFGenerator(xs, S1, S1))

    for dof in dg0.generate():
        assert dof.eval(lambda: 1) == 1

    converter = ElementSerialiser()
    encoded = converter.encode(dg0)
    decoded = converter.decode(encoded)

    converter = ElementSerialiser()
    encoded = converter.encode(edge)
    decoded = converter.decode(encoded)

    xs = [DOF(DeltaPairing(), PointKernel((-1,)))]
    dg1 = ElementTriple(decoded, (P1, CellL2, C0), DOFGenerator(xs, S2, S1))

    for dof in dg1.generate():
        assert np.allclose(abs(dof.eval(lambda x: x)), 1)


def test_repeated_objs():
    repeated_edge = Point(1, [vert, vert], vertex_num=2)
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
