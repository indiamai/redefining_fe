# from firedrake import *
from groups.groups import Group, S1, S2, S3
from cell_complex.cells import Point, Edge
import sympy as sp
from FIAT.functional import PointEvaluation
from FIAT.reference_element import UFCInterval, UFCTriangle
from spaces.continuous_spaces import ContinuousFunctionSpace


class Triple():
    
    def __init__(self, c, v, e):
        assert isinstance(c, Point)
        assert isinstance(e, E)
        # for space in v:
        #     assert isinstance(space, ContinuousFunctionSpace)
        self.C = c
        self.V = v
        self.E = e

    def generate(self):
        L = self.E.generate()
        return L

    def __iter__(self):
        yield self.C
        yield self.V
        yield self.E


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
                generated = l_g(g)
                if not isinstance(generated, list):
                    generated = [generated]
                ls.extend(generated)
        return ls


def immerse(g, G, triple):
    assert isinstance(triple, Triple)
    C, V, E = triple
    new_dofs = []
    for l in E.generate():
        new_dofs.append(trace(lambda x: g(G(x)), V, l))
    return new_dofs


def trace(G, V, v):
    # limited to point evaluations for now
    l_pts = v.pt_dict
    v_tilde_res = PointEvaluation(v.ref_el, G(list(l_pts.keys())[0]))
    return V[1].pullback(v_tilde_res)
