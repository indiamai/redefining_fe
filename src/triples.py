# from firedrake import *
from groups.groups import Group, S1, S2, S3
from cell_complex.cells import Point, Edge
import sympy as sp
from FIAT.functional import PointEvaluation
from FIAT.reference_element import UFCInterval, UFCTriangle
from spaces.element_sobolev_spaces import ElementSobolevSpace


class ElementTriple():

    def __init__(self, cell, spaces, dof_gen):
        assert isinstance(cell, Point)
        assert isinstance(dof_gen, DOFGenerator)

        self.cell = cell
        self.spaces = spaces
        self.DOFGenerator = dof_gen

    def generate(self):
        return self.DOFGenerator.generate()

    def __iter__(self):
        yield self.cell
        yield self.spaces
        yield self.DOFGenerator
    


class DOFGenerator():

    def __init__(self, generator_funcs, gen_group, trans_group):
        # assert isinstance(G_1, Group)
        # assert isinstance(G_2, Group)
        self.x = generator_funcs
        self.g1 = gen_group
        self.g2 = trans_group

    def __iter__(self):
        yield self.x
        yield self.g1
        yield self.g2

    def generate(self):
        ls = []
        for g in self.g1.members:
            for l_g in self.x:
                generated = l_g(g)
                if not isinstance(generated, list):
                    generated = [generated]
                ls.extend(generated)
        return ls


def immerse(g, attachment, triple):
    assert isinstance(triple, ElementTriple)

    C, V, DOFGenerator = triple
    new_dofs = []
    for generated_func in DOFGenerator.generate():
        new_dofs.append(trace(lambda x: g(attachment(x)), V, generated_func))
    return new_dofs


def trace(attachment, V, v):
    # should this be a property of the space
    # limited to point evaluations for now
    l_pts = v.pt_dict
    v_tilde_res = PointEvaluation(v.ref_el, attachment(list(l_pts.keys())[0]))
    return V[1].pullback(v_tilde_res)
