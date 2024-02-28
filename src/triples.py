# from firedrake import *
from groups.groups import Group, S1, S2, S3
from cell_complex.cells import Point, Edge
import sympy as sp
from FIAT.functional import PointEvaluation
from FIAT.reference_element import UFCInterval, UFCTriangle
from spaces.element_sobolev_spaces import ElementSobolevSpace
import matplotlib.pyplot as plt


class ElementTriple():

    def __init__(self, cell, spaces, dof_gen):
        assert isinstance(cell, Point)
        if isinstance(dof_gen, DOFGenerator):
            dof_gen = [dof_gen]
        for d in dof_gen:
            assert isinstance(d, DOFGenerator)

        self.cell = cell
        self.spaces = spaces
        self.DOFGenerator = dof_gen

        assert self.num_dofs() > self.spaces[0].subdegree

    def generate(self):
        res = []
        for dof_gen in self.DOFGenerator:
            res.extend(dof_gen.generate())
        return res

    def __iter__(self):
        yield self.cell
        yield self.spaces
        yield self.DOFGenerator

    def num_dofs(self):
        return sum([dof_gen.num_dofs() for dof_gen in self.DOFGenerator])

    def plot(self):
        dofs = self.generate()
        self.cell.plot(show=False, plain=True)
        for dof in dofs:
            l_pts = dof.pt_dict
            coord = list(l_pts.keys())[0]
            if self.cell.dimension == 1:
                coord = (coord, 0)
            if self.cell.dimension == 0:
                coord = (0, 0)
            plt.scatter(coord[0], coord[1], marker="o", color="blue")
        plt.show()


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

    def num_dofs(self):
        return len(self.x) * self.g1.size()

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

    C, V, _ = triple
    new_dofs = []
    for generated_func in triple.generate():
        new_dofs.append(trace(lambda x: g(attachment(x)), V, generated_func))
    return new_dofs


def trace(attachment, V, v):
    # should this be a property of the space
    # limited to point evaluations for now
    # how attachments work in 3d will require thought
    l_pts = v.pt_dict
    v_tilde_res = PointEvaluation(v.ref_el, attachment(list(l_pts.keys())[0]))
    return V[1].pullback(v_tilde_res)
