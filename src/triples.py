# from firedrake import *
from groups.groups import Group, S1, S2, S3
from cell_complex.cells import Point, Edge
import sympy as sp


class Triple():
    
    def __init__(self, c, v, e):
        assert isinstance(c, Point)
        assert isinstance(e, E)
        # for space in v:
        #     assert isinstance(space, FunctionSpace)
        self.C = c
        self.V = v
        self.E = e

    def generate(self):
        L = self.E.generate()
        return L


class E():

    def __init__(self, X, G_1, G_2):
        print(G_1)
        assert isinstance(G_1, Group)
        assert isinstance(G_2, Group)
        self.g1 = G_1
        self.g2 = G_2
        self.x = X

    def generate(self):
        l = []
        for g in self.g1.members:
            for l_g in self.x:
                l.append(l_g(g))
        return l
