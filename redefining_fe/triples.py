from redefining_fe.cells import Point
from redefining_fe.spaces.element_sobolev_spaces import ElementSobolevSpace, CellHCurl, CellHDiv
from redefining_fe.dof import DeltaPairing, L2InnerProd, MyTestFunction, PointKernel
from redefining_fe.traces import Trace
# from redefining_fe.serialisation import FETripleEncoder
from FIAT.dual_set import DualSet
from FIAT.finite_element import CiarletElement
import matplotlib.pyplot as plt
import inspect
import math
import finat.ufl
import warnings
import json


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
        id_counter = 0
        for dof_gen in self.DOFGenerator:
            generated = dof_gen.generate(self.cell, self.spaces[1], id_counter)
            res.extend(generated)
            id_counter += len(generated)
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

    def get_value_shape(self):
        # TODO Shape should be specificed somewhere else probably
        if self.spaces[0].vec:
            return (self.cell.get_spatial_dimension(),)
        else:
            return ()

    def to_ufl_elem(self):
        return IndiaTripleUFL(self)

    def to_fiat_elem(self):
        ref_el = self.cell.to_fiat()
        dofs = self.generate()
        degree = self.spaces[0].degree()
        entity_ids = {}
        entity_perms = {}
        nodes = []
        top = ref_el.get_topology()
        min_ids = self.cell.get_starter_ids()

        for dim in sorted(top):
            entity_ids[dim] = {i: [] for i in top[dim]}
            entity_perms[dim] = {}
            # perms = {0: [0]} if dim == 0 else self.make_entity_permutations(dim, degree - dim)
            # for entity in sorted(top[dim]):
            #         entity_perms[dim][entity] = perms

        entity_perms = None
        # print("DOFs", len(dofs))
        for i in range(len(dofs)):
            entity = dofs[i].trace_entity
            dim = entity.dim()
            entity_ids[dim][entity.id - min_ids[dim]].append(i)
            # print(dofs[i].id)
            nodes.append(dofs[i].convert_to_fiat(ref_el))
            # print(nodes[i].pt_dict)
        print("my ent ids", entity_ids)
        entity_perms = self.make_entity_permutations(self.cell.dim(), entity_ids, min_ids, dofs)

        form_degree = 1 if self.spaces[0].vec else 0
        dual = DualSet(nodes, ref_el, entity_ids, entity_perms)
        poly_set = self.spaces[0].to_ON_polynomial_set(ref_el)
        return CiarletElement(poly_set, dual, degree, form_degree)

    def make_entity_permutations_old(self, dim, entity_ids, min_ids):
        # limited to point eval
        # TODO: make this do the right thing
        # if npoints <= 0:
        #     return {o: [] for o in range(math.factorial(dim + 1))}
        id_counter = 0

        dof_info = {i: {"group": None, "ids": [], "dim": 0} for i in range(len(self.DOFGenerator))}
        for i in range(len(self.DOFGenerator)):
            sub_dofs = self.DOFGenerator[i].generate(self.cell, self.spaces[1], id_counter)
            dof_info[i]["dims"] = [sub_dofs[j].trace_entity.dim() for j in range(len(sub_dofs))]
            dof_info[i]["ents"] = [sub_dofs[j].trace_entity.id - min_ids[sub_dofs[j].trace_entity.dim()] for j in range(len(sub_dofs))]
            dof_info[i]["group"] = [sub_dofs[j].generation[dof_info[i]["dims"][j]].group for j in range(len(sub_dofs))]
            dof_info[i]["ids"] = list(range(id_counter, id_counter + len(sub_dofs)))
            id_counter += len(sub_dofs)

        res = {d: {} for d in range(dim + 1)}
        for d in range(dim + 1):
            for ent in entity_ids[d].keys():
                ent_dofs = entity_ids[d][ent]
                # min_ent = min(ent_dofs)
                # ent_dofs = [e - min_ent for e in ent_dofs]
                res[d][ent] = {o: ent_dofs[:] for o in range(math.factorial(d + 1))}
        print(dof_info)
        for i in range(len(self.DOFGenerator)):
            for j in range(len(dof_info[i]["ids"])):
                # base_val=min(dof_info[i]["ids"])
                orientation_reps = dof_info[i]["group"][j].compute_num_reps()
                dof_dim = dof_info[i]["dims"][j]
                dof_ent = dof_info[i]["ents"][j]
                dof_id = dof_info[i]["ids"][j]
                index = entity_ids[dof_dim][dof_ent].index(dof_id)
                # print(res)
                # print("orps", orientation_reps)
                for o in orientation_reps.keys():
                    res[dof_dim][dof_ent][o][index] = orientation_reps[o]
        print(res)
        return res
        # raise NotImplementedError("TODO work out orientations")

    def make_entity_permutations(self, dim, entity_ids, min_ids, dofs):
        res = {d: {} for d in range(dim + 1)}
        for d in range(dim + 1):
            for ent in entity_ids[d].keys():
                ent_dofs = entity_ids[d][ent]
                if len(ent_dofs) > 0:
                    min_ent = min(ent_dofs)
                    ent_dofs = [e - min_ent for e in ent_dofs]
                res[d][ent] = {o: ent_dofs[:] for o in range(math.factorial(d + 1))}
        for dof in dofs:
            dof_dim = dof.trace_entity.dim()
            dof_ent = dof.trace_entity.id - min_ids[dof_dim]
            group = dof.generation[dof_dim].group
            dof_index = entity_ids[dof_dim][dof_ent].index(dof.id)
            dof_group = None
            for gen in self.DOFGenerator:
                if dof.id in gen.dof_ids:
                    dof_group = gen.dof_ids
            assert dof_group is not None
            orientation_reps = group.compute_num_reps(dof_group)
            for o in orientation_reps.keys():
                res[dof_dim][dof_ent][o][dof_index] = orientation_reps[o][dof.sub_id]
            # need to maintain structure of group generation \
            # - what are the overall ids of the dofs generated by each group
        print(res)
        return res

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
                ax.text(*coord, dof.id)

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
                ax.text(*coord, dof.id)

            plt.show()
        else:
            raise ValueError("Plotting not supported in this dimension")
        
    # def to_json(self, filename="triple.json"):
        # encoded = json.dumps(self, cls=FETripleEncoder)
        # print(encoded)
        # with open(filename, "w+") as f:
        #     f.write(encoded)





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

    def generate(self, cell, space, id_counter):
        if self.ls is None:
            self.ls = []
            for l_g in self.x:
                i = 0
                for g in self.g1.members():
                    generated = l_g(g)
                    if not isinstance(generated, list):
                        generated = [generated]
                    for dof in generated:
                        dof.add_context(cell, space, g, id_counter, i)
                        id_counter += 1
                        i += 1
                    self.ls.extend(generated)
            self.dof_numbers = len(self.ls)
            self.dof_ids = [dof.id for dof in self.ls]
            return self.ls
        return self.ls

    def generate_by_x(self, cell, space, id_counter):
        res_ls = []
        for l_g in self.x:
            ls = []
            for g in self.g1.members():
                generated = l_g(g)
                if not isinstance(generated, list):
                    generated = [generated]
                for dof in generated:
                    dof.add_context(cell, space, g, id=id_counter)
                    id_counter += 1
                ls.extend(generated)
            res_ls.append(ls)

        return res_ls

    def __repr__(self):
        repr_str = ""
        for x_elem in self.x:
            repr_str += "g(" + str(x_elem) + ")"
        return repr_str



class ImmersedDOFs():

    def __init__(self, target_cell, triple, trace, start_node=0):
        self.target_cell = target_cell
        self.triple = triple
        self.C, self.V, self.E = triple
        self.trace = trace(target_cell)
        self.start_node = start_node

    def __call__(self, g):
        target_node, o = self.target_cell.permute_entities(g, self.C.dim())[self.start_node]
        if self.C.dim() > 0 and o != o.group.identity:
            return []
        attachment = self.target_cell.cell_attachment(target_node)
        new_dofs = []

        def oriented_attachment(*x):
            return attachment(*o(x))

        for generated_dof in self.triple.generate():
            new_dof = generated_dof.immerse(self.target_cell.get_node(target_node),
                                            oriented_attachment,
                                            self.trace, g, self.triple)
            new_dofs.append(new_dof)
        return new_dofs

    def __repr__(self):
        repr_str = ""
        for dof_gen in self.E:
            repr_str += "Im_" + str(self.trace) + "_" + str(self.target_cell) + "(" + str(dof_gen) + ")"
        return repr_str


def immerse(target_cell, triple, target_space, node=0):
    return ImmersedDOFs(target_cell, triple, target_space, node)


class IndiaTripleUFL(finat.ufl.FiniteElementBase):
    """
    TODO: Need to deal with cases where value shape and reference value shape are different
    """

    def __init__(self, triple, cell=None):
        self.triple = triple
        if not cell:
            cell = self.triple.cell.to_ufl()

        # this isn't really correct
        degree = self.triple.spaces[0].degree()

        super(IndiaTripleUFL, self).__init__("IT", cell, degree, None, triple.get_value_shape(), triple.get_value_shape())

    def __repr__(self):
        return "FiniteElement(%s, (%s, %s, %s), %s)" % (
            repr(self.triple.cell), repr(self.triple.spaces[0]), repr(self.triple.spaces[1]), repr(self.triple.spaces[2]), "X")

    def __str__(self):
        return "<Custom%sElem on %s>" % (self.triple.spaces[0], self.triple.cell)

    def mapping(self):
        if isinstance(self.sobolev_space, CellHCurl):
            return "covariant Piola"
        elif isinstance(self.sobolev_space, CellHDiv):
            return "contravariant Piola"
        else:
            return "identity"

    def sobolev_space(self):
        return self.triple.spaces[1]

    def reconstruct(self, family=None, cell=None, degree=None, quad_scheme=None, variant=None):
        warnings.warn("Modifying FE triple")
        return IndiaTripleUFL(self.triple, cell=cell)
