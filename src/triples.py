# from firedrake import *
from groups.groups import Group, S1, S2, S3
from cell_complex.cells import Point, Edge
import sympy as sp
from FIAT.functional import PointEvaluation
from FIAT.reference_element import UFCInterval, UFCTriangle


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
    
    def split(self):
        return self.C, self.V, self.E


class E():

    def __init__(self, X, G_1, G_2):
        assert isinstance(G_1, Group)
        assert isinstance(G_2, Group)
        self.g1 = G_1
        self.g2 = G_2
        self.x = X

    def generate(self):
        ls = []
        for g in self.g1.members:
            for l_g in self.x:
                ls.append(l_g(g))
        return ls


def immerse(g, G, triple):
    # this is limited to point evaluations for now.
    # missing trace
    assert isinstance(triple, Triple)
    C, V, E = triple.split()
    new_dofs = []
    for l in E.generate():
        l_pts = l.pt_dict
        print(list(l_pts.keys())[0])
        new_dofs.append(PointEvaluation(l.ref_el, G(g(list(l_pts.keys())[0]))))
    print("generated")
    print(new_dofs)
    for n in new_dofs:
        print(n.tostr())
    return new_dofs
