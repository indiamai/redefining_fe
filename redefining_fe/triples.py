from redefining_fe.cells import Point
from redefining_fe.spaces.element_sobolev_spaces import ElementSobolevSpace, CellHCurl, CellHDiv
from redefining_fe.dof import DeltaPairing, L2InnerProd, MyTestFunction, PointKernel
from redefining_fe.traces import Trace
from redefining_fe.groups import perm_matrix_to_perm_array
from FIAT.dual_set import DualSet
from FIAT.finite_element import CiarletElement
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import inspect
import finat.ufl
import warnings
import numpy as np
import scipy


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
        if self.spaces[0].set_shape:
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
        poly_set = self.spaces[0].to_ON_polynomial_set(ref_el)

        for dim in sorted(top):
            entity_ids[dim] = {i: [] for i in top[dim]}
            entity_perms[dim] = {}

        for i in range(len(dofs)):
            entity = dofs[i].trace_entity
            dim = entity.dim()
            entity_ids[dim][entity.id - min_ids[dim]].append(i)
            nodes.append(dofs[i].convert_to_fiat(ref_el, degree))
        entity_perms, pure_perm = self.make_dof_perms(ref_el, entity_ids, nodes, poly_set)

        form_degree = 1 if self.spaces[0].set_shape else 0
        print("my", [n.pt_dict for n in nodes])
        # TODO: Change this when Dense case in Firedrake
        if pure_perm:
            dual = DualSet(nodes, ref_el, entity_ids, entity_perms)
        else:
            dual = DualSet(nodes, ref_el, entity_ids)
        return CiarletElement(poly_set, dual, degree, form_degree)

    def plot(self, filename="temp.png"):
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

            fig.savefig(filename)
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

            fig.savefig(filename)
        else:
            raise ValueError("Plotting not supported in this dimension")

    def compute_dense_matrix(self, ref_el, entity_ids, nodes, poly_set):
        dual = DualSet(nodes, ref_el, entity_ids)

        old_coeffs = poly_set.get_coeffs()
        dualmat = dual.to_riesz(poly_set)

        shp = dualmat.shape
        A = dualmat.reshape((shp[0], -1))
        B = old_coeffs.reshape((shp[0], -1))
        V = np.dot(A, np.transpose(B))

        with warnings.catch_warnings():
            warnings.filterwarnings("error")
            try:
                new_coeffs_flat = scipy.linalg.solve(V, B, transposed=True)
            except (scipy.linalg.LinAlgWarning, scipy.linalg.LinAlgError):
                raise np.linalg.LinAlgError("Singular Vandermonde matrix")
        return A, new_coeffs_flat

    def make_overall_dense_matrices(self, ref_el, entity_ids, nodes, poly_set):
        min_ids = self.cell.get_starter_ids()
        dim = self.cell.dim()
        e = self.cell
        e_id = e.id - min_ids[dim]
        res_dict = {dim: {e_id: {}}}
        degree = self.spaces[0].degree()
        original_V, original_basis = self.compute_dense_matrix(ref_el, entity_ids, nodes, poly_set)
        for g in self.cell.group.members():
            val = g.numeric_rep()
            if g.perm.is_Identity:
                res_dict[dim][e_id][val] = np.eye(len(nodes))
            else:
                new_nodes = [d(g).convert_to_fiat(ref_el, degree) for d in self.generate()]
                transformed_V, transformed_basis = self.compute_dense_matrix(ref_el, entity_ids, new_nodes, poly_set)
                res_dict[dim][e_id][val] = np.matmul(transformed_basis, original_V.T)
        return res_dict

    def make_dof_perms(self, ref_el, entity_ids, nodes, poly_set):
        dofs = self.generate()
        min_ids = self.cell.get_starter_ids()
        entity_associations = {dim: {e.id - min_ids[dim]: {} for e in self.cell.d_entities(dim, get_class=True)}
                               for dim in range(self.cell.dim() + 1)}
        cell_dim = self.cell.dim()
        cell_dict = entity_associations[cell_dim][self.cell.id - min_ids[cell_dim]]
        pure_perm = True

        # construct mapping of entities to the dof generators and the dofs they generate
        for d in dofs:
            sub_dim = d.trace_entity.dim()
            sub_dict = entity_associations[sub_dim][d.trace_entity.id - min_ids[sub_dim]]
            dof_gen = str(d.generation[sub_dim])
            if not len(d.generation[sub_dim].g2.members()) == 1:
                pure_perm = False
            if dof_gen in sub_dict.keys():
                sub_dict[dof_gen] += [d]
            else:
                sub_dict[dof_gen] = [d]
            if sub_dim != cell_dim:
                dof_gen = str(d.generation[cell_dim])
                if not len(d.generation[cell_dim].g2.members()) == 1:
                    pure_perm = False

                if dof_gen in cell_dict.keys():
                    cell_dict[dof_gen] += [d]
                else:
                    cell_dict[dof_gen] = [d]

        if pure_perm is False:
            # TODO think about where this call goes
            return self.make_overall_dense_matrices(ref_el, entity_ids, nodes, poly_set), pure_perm

        dof_id_mat = np.eye(len(dofs))
        oriented_mats_by_entity = {}
        flat_by_entity = {}

        # for each entity, look up generation on that entity and permute the
        # dof mapping according to the generation
        for dim in range(self.cell.dim()):
            oriented_mats_by_entity[dim] = {}
            flat_by_entity[dim] = {}
            ents = self.cell.d_entities(dim, get_class=True)
            for e in ents:
                e_id = e.id - min_ids[dim]
                members = e.group.members()
                oriented_mats_by_entity[dim][e_id] = {}
                flat_by_entity[dim][e_id] = {}
                for g in members:
                    val = g.numeric_rep()
                    oriented_mats_by_entity[dim][e_id][val] = dof_id_mat.copy()
                    flat_by_entity[dim][e_id][val] = []
                    for dof_gen in entity_associations[dim][e_id].keys():
                        ent_dofs = entity_associations[dim][e_id][dof_gen]
                        ent_dofs_ids = np.array([ed.id for ed in ent_dofs], dtype=int)
                        dof_gen_class = ent_dofs[0].generation[dim]
                        # TODO not sure about correctness of this
                        if g.perm.is_Identity or (pure_perm and len(ent_dofs_ids) == 1):
                            oriented_mats_by_entity[dim][e_id][val][np.ix_(ent_dofs_ids, ent_dofs_ids)] = np.eye(len(ent_dofs_ids))
                            flat_by_entity[dim][e_id][val] = perm_matrix_to_perm_array(np.eye(len(ent_dofs_ids)))
                        elif g in dof_gen_class.g1.members():
                            sub_mat = g.matrix_form()
                            oriented_mats_by_entity[dim][e_id][val][np.ix_(ent_dofs_ids, ent_dofs_ids)] = sub_mat.copy()
                            flat_by_entity[dim][e_id][val] = perm_matrix_to_perm_array(sub_mat)

        oriented_mats_overall = {}
        dim = self.cell.dim()
        e = self.cell
        e_id = e.id - min_ids[dim]
        flat_by_entity[dim] = {e_id: {}}
        members = e.group.members()
        for g in members:
            val = g.numeric_rep()
            oriented_mats_overall[val] = dof_id_mat.copy()
            flat_by_entity[dim][e_id][val] = []
            for dof_gen in entity_associations[dim][e_id].keys():
                ent_dofs = entity_associations[dim][e_id][dof_gen]
                ent_dofs_ids = np.array([ed.id for ed in ent_dofs], dtype=int)
                dof_gen_class = ent_dofs[0].generation
                for key in dof_gen_class.keys():
                    if not key == dim:
                        immersed_dim = key
                        for sub_e, sub_g in e.permute_entities(g, immersed_dim):
                            sub_e = e.get_node(sub_e)
                            sub_e_id = sub_e.id - min_ids[sub_e.dim()]
                            sub_ent_assoc = entity_associations[immersed_dim][sub_e_id][str(dof_gen_class[immersed_dim])]
                            sub_mat = oriented_mats_by_entity[immersed_dim][sub_e_id][sub_g.numeric_rep()][np.ix_(ent_dofs_ids, ent_dofs_ids)]

                            g_sub_mat = g.matrix_form()
                            expanded = np.kron(g_sub_mat, np.eye(len(sub_ent_assoc)))
                            oriented_mats_overall[val][np.ix_(ent_dofs_ids, ent_dofs_ids)] = np.matmul(sub_mat, expanded).copy()
                    elif len(dof_gen_class.keys()) == 1:
                        if g in dof_gen_class[dim].g1.members():
                            sub_mat = g.matrix_form()
                            oriented_mats_overall[val][np.ix_(ent_dofs_ids, ent_dofs_ids)] = sub_mat.copy()
                        elif pure_perm and len(dof_gen_class[dim].g1.members()) > 1:
                            sub_mat = g.matrix_form()
                            oriented_mats_overall[val][np.ix_(ent_dofs_ids, ent_dofs_ids)] = sub_mat.copy()

        for val, mat in oriented_mats_overall.items():
            cell_dofs = entity_ids[dim][0]
            flat_by_entity[dim][e_id][val] = perm_matrix_to_perm_array(mat[np.ix_(cell_dofs, cell_dofs)])

        return flat_by_entity, pure_perm

    def _to_dict(self):
        o_dict = {"cell": self.cell, "spaces": self.spaces, "dofs": self.DOFGenerator}
        return o_dict

    def dict_id(self):
        return "Triple"

    def _from_dict(o_dict):
        return ElementTriple(o_dict["cell"], o_dict["spaces"], o_dict["dofs"])


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
                        dof.add_context(self, cell, space, g, id_counter, i)
                        id_counter += 1
                        i += 1
                    self.ls.extend(generated)
            self.dof_numbers = len(self.ls)
            self.dof_ids = [dof.id for dof in self.ls]
        return self.ls

    # def generate_by_x(self, cell, space, id_counter):
    #     res_ls = []
    #     for l_g in self.x:
    #         ls = []
    #         for g in self.g1.members():
    #             generated = l_g(g)
    #             if not isinstance(generated, list):
    #                 generated = [generated]
    #             for dof in generated:
    #                 dof.add_context(self, cell, space, g, id=id_counter)
    #                 id_counter += 1
    #             ls.extend(generated)
    #         res_ls.append(ls)

        # return res_ls

    def __repr__(self):
        repr_str = "DOFGen("
        for x_elem in self.x:
            repr_str += "g(" + str(x_elem) + ")"
        repr_str += str(self.g1) + str(self.g2) + ")"
        return repr_str

    def _to_dict(self):
        o_dict = {"x": self.x, "g1": self.g1, "g2": self.g2}
        return o_dict

    def dict_id(self):
        return "DOFGen"

    def _from_dict(obj_dict):
        return DOFGenerator(obj_dict["x"], obj_dict["g1"], obj_dict["g2"])


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
            raise ValueError("Not matching orientation - groups incorrect")
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

    def _to_dict(self):
        o_dict = {"target_cell": self.target_cell, "triple": self.triple, "trace": self.trace}
        return o_dict

    def dict_id(self):
        return "ImmersedDOF"

    def _from_dict(obj_dict):
        return ImmersedDOFs(obj_dict["target_cell"], obj_dict["triple"], obj_dict["trace"])


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

        super(IndiaTripleUFL, self).__init__("IT", cell, degree, None, triple.get_value_shape())

    def __repr__(self):
        return "FiniteElement(%s, %s, (%s, %s, %s), %s)" % (
               repr(self.triple.DOFGenerator), repr(self.triple.cell), repr(self.triple.spaces[0]), repr(self.triple.spaces[1]), repr(self.triple.spaces[2]), "X")

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
