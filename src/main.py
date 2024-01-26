# top level file for linking together all the packages
from firedrake import *
from groups.group import S1, S2, S3
from cell_complex.cells import *
e = S3()

print(e)

edge = Interval()
edge_cell = edge.default_cell_complex()


def generate(C, V, E):
    (dofs, g1, g2) = E
    dofs_under_g1 = []
    for g in g1:
        dofs_under_g1.extend([dof(g) for dof in dofs])
    # then replace all immersion operators with the right ones.
