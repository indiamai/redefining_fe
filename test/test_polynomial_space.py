from redefining_fe import *
import sympy as sp
from redefining_fe.spaces.polynomial_spaces import PolynomialSpace, RestrictedPolynomialSpace, ConstructedPolynomialSpace
from FIAT.polynomial_set import ONPolynomialSet


def test_instantiation():
    cell = Point(0)
    x = sp.Symbol("x")
    composite = P0 + x*P0

    assert isinstance(composite.to_ON_polynomial_set(cell), ONPolynomialSet)


def test_restriction():
    restricted = P3.restrict(2, 3)
    # doesn't contain constants
    assert restricted.subdegree == -1
    assert restricted.superdegree == 3

    x = sp.Symbol("x")
    composite = P2 + x*restricted
    assert isinstance(composite, ConstructedPolynomialSpace)
    print(composite.spaces)

    not_restricted = P3.restrict(0, 3)
    assert isinstance(not_restricted, PolynomialSpace)
    assert not isinstance(not_restricted, RestrictedPolynomialSpace)
