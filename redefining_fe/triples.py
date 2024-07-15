from redefining_fe.cells import Point
from redefining_fe.spaces.element_sobolev_spaces import ElementSobolevSpace
from redefining_fe.dof import DeltaPairing, L2InnerProd, MyTestFunction, PointKernel
from redefining_fe.traces import Trace
import matplotlib.pyplot as plt
import inspect


class ElementTriple():
    """
    Class to represent the three core parts of the element

    :param: cell: CellComplex
    :param: spaces: Triple of spaces: (PolynomialSpace, SobolovSpace, InterpolationSpace)
    :param: dof_gen: Generator Triple to generate the degrees of freedom.
    """

    def __init__(self, cell, spaces, dof_gen):
        assert isinstance(cell, Point)
        if isinstance(dof_gen, DOFGenerator):
            dof_gen = [dof_gen]
        for d in dof_gen:
            assert isinstance(d, DOFGenerator)
            d.add_cell(cell)

        self.cell = cell
        cell_spaces = []
        for space in spaces:
            # TODO: Fix this to a more sensible condition when all spaces
            # implemented
            if inspect.isclass(space) and issubclass(space, ElementSobolevSpace):
                cell_spaces.append(space(cell))
            else:
                cell_spaces.append(space)
        self.spaces = tuple(cell_spaces)
        self.DOFGenerator = dof_gen

    def generate(self):
        res = []
        for dof_gen in self.DOFGenerator:
            res.extend(dof_gen.generate(self.cell, self.spaces[1]))
        return res

    def __iter__(self):
        yield self.cell
        yield self.spaces
        yield self.DOFGenerator

    def num_dofs(self):
        return sum([dof_gen.num_dofs() for dof_gen in self.DOFGenerator])

    def get_dof_info(self, dof):
        if dof.trace_entity.dimension == 0:
            center = self.cell.cell_attachment(dof.trace_entity.id)()
            color = "b"
        elif dof.trace_entity.dimension == 1:
            color = "r"
            center = self.cell.cell_attachment(dof.trace_entity.id)(0)
        elif dof.trace_entity.dimension == 2:
            color = "g"
            center = self.cell.cell_attachment(dof.trace_entity.id)(0, 0)
        else:
            color = "b"
            center = None

        return center, color

    def plot(self):
        # point evaluation nodes only
        dofs = self.generate()
        identity = MyTestFunction(lambda *x: x)

        if self.cell.dimension == 0:
            raise ValueError(" Dimension 0 cells cannot be plotted")

        if self.cell.dimension < 3:
            fig = plt.figure()
            ax = plt.gca()
            self.cell.plot(show=False, plain=True, ax=ax)
            for dof in dofs:
                center, color = self.get_dof_info(dof)
                if isinstance(dof.pairing, DeltaPairing) and isinstance(dof.kernel, PointKernel):
                    coord = dof.eval(identity, pullback=False)
                elif isinstance(dof.pairing, L2InnerProd):
                    coord = center
                if len(coord) == 1:
                    coord = (coord[0], 0)
                if isinstance(dof.target_space, Trace):
                    dof.target_space.plot(ax, coord, dof.trace_entity, dof.g, color=color)
                else:
                    ax.scatter(*coord, color=color)

            plt.show()
        elif self.cell.dimension == 3:
            fig = plt.figure()
            ax = fig.add_subplot(projection='3d')
            self.cell.plot3d(show=False, ax=ax)
            for dof in dofs:
                center, color = self.get_dof_info(dof)
                if center is None:
                    center = [0, 0, 0]
                if isinstance(dof.pairing, DeltaPairing):
                    coord = dof.eval(identity, pullback=False)
                    if isinstance(dof.target_space, Trace):
                        dof.target_space.plot(ax, coord, dof.trace_entity, dof.g, color=color)
                    else:
                        ax.scatter(*coord, color=color)
                elif isinstance(dof.pairing, L2InnerProd):
                    dof.target_space.plot(ax, center, dof.trace_entity, dof.g, color=color, length=0.2)

            plt.show()
        else:
            raise ValueError("Plotting not supported in this dimension")


class DOFGenerator():

    def __init__(self, generator_funcs, gen_group, trans_group):
        # assert isinstance(G_1, Group)
        # assert isinstance(G_2, Group)
        self.x = generator_funcs
        self.g1 = gen_group
        self.g2 = trans_group
        self.dof_numbers = None
        self.ls = None

    def __iter__(self):
        yield self.x
        yield self.g1
        yield self.g2

    def add_cell(self, cell):
        self.g1 = self.g1.add_cell(cell)
        self.g2 = self.g2.add_cell(cell)

    def num_dofs(self):
        if self.dof_numbers is None:
            raise ValueError("DOFs not generated yet")
        return self.dof_numbers

    def generate(self, cell, space):
        if self.ls is None:
            self.ls = []
            for g in self.g1.members():
                for l_g in self.x:
                    generated = l_g(g)
                    if not isinstance(generated, list):
                        generated = [generated]
                    for dof in generated:
                        dof.add_context(cell, space)
                    self.ls.extend(generated)
            self.dof_numbers = len(self.ls)
            return self.ls
        return self.ls

    def __repr__(self):
        repr_str = ""
        for x_elem in self.x:
            repr_str += "g(" + str(x_elem) + ")"
        return repr_str


class ImmersedDOF():

    def __init__(self, target_cell, triple, target_space, start_node=0):
        self.target_cell = target_cell
        self.triple = triple
        self.C, self.V, self.E = triple
        self.target_space = target_space(target_cell)
        self.start_node = start_node

    def __call__(self, g):
        # print(self.start_node)
        # print(g)
        # print(self.target_cell.permute_entities(g, self.C.dim()))
        target_node, o = self.target_cell.permute_entities(g, self.C.dim())[self.start_node]

        # print("Attaching node", target_node)
        # print("Orientation", o)
        if self.C.dim() > 0 and o != o.group.identity:
            # print("orientation doesn't match")
            return []
        attachment = self.target_cell.cell_attachment(target_node)
        new_dofs = []

        def oriented_attachment(*x):
            return attachment(*o(x))

        for generated_dof in self.triple.generate():
            new_dof = generated_dof.immerse(self.target_cell.get_node(target_node),
                                            oriented_attachment,
                                            self.target_space, g)
            new_dofs.append(new_dof)
        return new_dofs

    def __repr__(self):
        repr_str = ""
        for dof_gen in self.E:
            repr_str += "Im_" + str(self.target_space) + "_" + str(self.target_cell) + "(" + str(dof_gen) + ")"
        return repr_str


def immerse(target_cell, triple, target_space, node=0):
    return ImmersedDOF(target_cell, triple, target_space, node)
