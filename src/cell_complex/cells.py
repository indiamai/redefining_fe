import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import ConvexHull
import itertools
import networkx as nx
import groups.new_groups
import copy
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d
from sympy.combinatorics.named_groups import SymmetricGroup, Permutation, PermutationGroup


class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        super().__init__((0, 0), (0, 0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def do_3d_projection(self, renderer=None):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, self.axes.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))

        return np.min(zs)


def symmetric_group_rep(d):
    return groups.new_groups.GroupRepresentation(SymmetricGroup(d))
    if d == 1:
        return groups.new_groups.S1
    if d == 2:
        return groups.new_groups.S2
    if d == 3: 
        return groups.new_groups.S3
    if d == 4:
        return groups.new_groups.S4


def topo_pos(G):
    """ helper function for hasse diagram visualisation
    Display in topological order, with simple offsetting for legibility"""
    pos_dict = {}
    for i, node_list in enumerate(nx.topological_generations(G)):
        x_offset = len(node_list) / 2
        y_offset = 0
        for j, name in enumerate(node_list):
            pos_dict[name] = (j - x_offset, i - j * y_offset)

    return pos_dict


def fold_reduce(func_list, *prev):
    """ nested function composition helper function, right to left """
    for func in reversed(func_list):
        # print(prev)
        prev = func(*prev)
    # print(prev)
    return prev


def normalise(v):
    norm = np.linalg.norm(v)
    return v / norm


def make_arrow(ax, mid, edge, direction=1):
    delta = 0.0001 if direction >= 0 else -0.0001
    x, y = edge(mid)
    dir_x, dir_y = edge(mid + delta)
    ax.arrow(x, y, dir_x-x, dir_y-y, head_width=0.05, head_length=0.1)


def make_arrow_3d(ax, mid, edge, direction=1):
    delta = 0.0001 if direction >= 0 else -0.0001
    x, y, z = edge(mid)
    dir_x, dir_y, dir_z = edge(mid + delta)
    a = Arrow3D([x, dir_x], [y, dir_y], [z, dir_z], mutation_scale=10, arrowstyle="-|>", color="black")
    ax.add_artist(a)


class Point():

    id_iter = itertools.count()

    def __init__(self, d, edges=[], oriented=False, group=None):
        self.id = next(self.id_iter)
        self.dimension = d
        if d == 0:
            assert (edges == [])

        self.oriented = oriented
        self.group = None
        self.G = nx.DiGraph()
        self.G.add_node(self.id, point_class=self)
        for edge in edges:
            assert edge.lower_dim() < self.dimension
            self.G.add_edge(self.id, edge.point.id, edge_class=edge)
        self.G = nx.compose_all([self.G]
                                + [edge.point.graph() for edge in edges])
        self.connections = edges
        group = self.compute_cell_group()
        if group:
            self.group = group.add_cell(self)

    def compute_cell_group(self):
        verts = self.vertices()
        v_coords = self.vertices(return_coords=True)
        coords_dict = {v: coord for v, coord in zip(verts, v_coords)}
        print(verts)
        n = len(verts)
        max_group = SymmetricGroup(n)
        edges = [edge.vertices() for edge in self.edges(get_class=True)]

        accepted_perms = max_group.elements
        print("max num", len(accepted_perms))
        if n > 2:
            for element in max_group.elements:
                reordered = element(verts)
                for edge in edges:
                    diff = np.subtract(v_coords[reordered.index(edge[0])], v_coords[reordered.index(edge[1])])
                    edge_len = np.sqrt(np.dot(diff, diff))
                    if not np.allclose(edge_len, 2):
                        accepted_perms.remove(element)
                        break
        print("Accepted", len(accepted_perms))
        print(accepted_perms)
        return groups.new_groups.GroupRepresentation(PermutationGroup(list(accepted_perms)))




    def dim(self):
        return self.dimension

    def graph_dim(self):
        if self.oriented:
            dim = self.dimension + 1
        else:
            dim = self.dimension
        return dim

    def graph(self):
        if self.oriented:
            temp_G = self.G.copy()
            temp_G.remove_node(-1)
            return temp_G
        return self.G

    def hasse_diagram(self, counter=0):
        ax = plt.axes()
        nx.draw_networkx(self.graph(), pos=topo_pos(self.graph()),
                         with_labels=True, ax=ax)
        plt.show()

    def d_entities(self, d, get_class=False):
        levels = [sorted(generation)
                  for generation in nx.topological_generations(self.G)]
        if get_class:
            return [self.G.nodes.data("point_class")[i]
                    for i in levels[self.graph_dim() - d]]
        return levels[self.graph_dim() - d]

    def get_node(self, node):
        return self.G.nodes.data("point_class")[node]

    def dim_of_node(self, node):
        levels = [sorted(generation)
                  for generation in nx.topological_generations(self.G)]
        for i in range(len(levels)):
            if node in levels[i]:
                return self.graph_dim() - i
        raise "Error: Node not found in graph"
    
    def vertices(self, get_class=False, return_coords=False):
        if self.oriented:
            verts = self.oriented.permute(self.d_entities(0, get_class))
        else:
            verts = self.d_entities(0, get_class)
        if return_coords:
            top_level_node = self.d_entities(self.graph_dim())[0]
            if self.dimension == 0:
                return [()]
            return [self.attachment(top_level_node, v)() for v in verts]
        return verts

    
    def edges(self, get_class=False):
        if self.oriented:
            return self.oriented.permute(self.d_entities(1, get_class))
        return self.d_entities(1, get_class)
    
    def permute_entities(self, g, d):
        verts = self.vertices()
        entities = self.d_entities(d)
        reordered = g.permute(verts)

        if d == 0:
            return list(zip(reordered, [lambda *x: () for r in reordered]))

        entity_dict = {}
        reordered_entity_dict = {}

        for ent in entities:
            entity_verts = []
            for v in verts:
                connected = list(nx.all_simple_edge_paths(self.G, ent, v))
                if len(connected) > 0:
                    entity_verts.append(v)
            entity_dict[ent] = tuple(entity_verts)
            reordered_entity_dict[ent] = tuple([reordered[i] for i in entity_verts])

        reordered_entities = [tuple() for e in range(len(entities))]
        min_id = min(entities)
        entity_vert_num = len(entity_dict[min_id])
        entity_group = symmetric_group_rep(entity_vert_num).add_cell(self.get_node(min_id))
        for ent in entities:
            for ent1 in entities:
                if set(entity_dict[ent]) == set(reordered_entity_dict[ent1]):
                    if entity_dict[ent] != reordered_entity_dict[ent1]:
                        o = entity_group.transform_between_perms(entity_dict[ent], reordered_entity_dict[ent1])
                        reordered_entities[ent1 - min_id] = (ent, o)
                    else:
                        reordered_entities[ent1 - min_id] = (ent, entity_group.identity)
        return reordered_entities


    def basis_vectors(self, return_coords=True, entity=None):
        if not entity:
            entity = self
        vertices = entity.vertices()
        if self.dim == 0:
            raise ValueError("Dimension 0 entities cannot have Basis Vectors")
        top_level_node = self.d_entities(self.graph_dim())[0]
        v_0 = vertices[0]
        if return_coords:
            v_0_coords = self.attachment(top_level_node, v_0)()
        basis_vecs = []
        for v in vertices[1:]:
            if return_coords:
                v_coords = self.attachment(top_level_node, v)()
                sub = normalise(np.subtract(v_coords, v_0_coords))
                if not hasattr(sub, "__iter__"):
                    basis_vecs.append((sub,))
                else:
                    basis_vecs.append(tuple(sub))
            else:
                basis_vecs.append((v, v_0))
        return basis_vecs


    def plot(self, show=True, plain=False, ax=None):
        """ for now into 2 dimensional space """

        top_level_node = self.d_entities(self.graph_dim())[0]
        xs = np.linspace(-1, 1, 20)
        if ax is None:
            ax = plt.gca()

        if self.dimension == 1:
            # line plot in 1D case
            nodes = self.d_entities(0)
            points = []
            for node in nodes:
                attach = self.attachment(top_level_node, node)
                points.extend(attach())
            plt.plot(np.array(points), np.zeros_like(points), color="black")


        for i in range(self.dimension - 1, -1, -1):
            nodes = self.d_entities(i)
            vert_coords = []
            for node in nodes:
                attach = self.attachment(top_level_node, node)
                if i == 0:
                    plotted = attach()
                    if len(plotted) < 2:
                        plotted = (plotted[0], 0)
                    vert_coords += [plotted]
                    if not plain:
                        plt.plot(plotted[0], plotted[1], 'bo')
                        plt.annotate(node, (plotted[0], plotted[1]))
                elif i == 1:
                    edgevals = np.array([attach(x) for x in xs])
                    if len(edgevals[0]) < 2:
                        plt.plot(edgevals[:, 0], 0, color="black")
                    else:
                        plt.plot(edgevals[:, 0], edgevals[:, 1], color="black")
                    if not plain:
                        make_arrow(ax, 0, attach)
                else:
                    raise ValueError("General plotting not implemented")
            # if i == 2:
            #     if len(vert_coords) > 2:
            #         hull = ConvexHull(vert_coords)
            #         plt.fill(vert_coords[hull.vertices, 0], vert_coords[hull.vertices, 1], alpha=0.5)
        if show:
            plt.show()

    def plot3d(self, show=True, ax = None):
        assert self.dimension == 3
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(projection='3d')
        xs = np.linspace(-1, 1, 20)

        top_level_node = self.d_entities(self.graph_dim())[0]
        nodes = self.d_entities(0)
        for node in nodes:
            attach = self.attachment(top_level_node, node)
            plotted = attach()
            ax.scatter(plotted[0], plotted[1], plotted[2],color="black")

        nodes = self.d_entities(1)
        for node in nodes:
            attach = self.attachment(top_level_node, node)
            edgevals = np.array([attach(x) for x in xs])
            ax.plot3D(edgevals[:, 0], edgevals[:, 1], edgevals[:, 2], color="black")
            make_arrow_3d(ax, 0, attach)
        if show:
            plt.show()

    def attachment(self, source, dst):
        if source == dst:
            return lambda *x: x
        
        paths = nx.all_simple_edge_paths(self.G, source, dst)
        attachments = [[self.G[s][d]["edge_class"]
                        for (s, d) in path] for path in paths]

        if len(attachments) == 0:
            raise ValueError("No paths from node {} to node {}"
                             .format(source, dst))
        # print(attachments)
        # check all attachments resolve to the same function
        if len(attachments) > 1:
            dst_dim = self.dim_of_node(dst)
            basis = np.eye(dst_dim)
            if dst_dim == 0:
                vals = [fold_reduce(attachment) for attachment in attachments]
                # for val in vals:
                #     print(val)
                assert all(np.isclose(val, vals[0]).all() for val in vals)
            else:
                for i in range(dst_dim):
                    vals = [fold_reduce(attachment, *tuple(basis[i].tolist()))
                            for attachment in attachments]
                    assert all(np.isclose(val, vals[0]).all() for val in vals)

        return lambda *x: fold_reduce(attachments[0], *x)

    def cell_attachment(self, dst):
        top_level_node = self.d_entities(self.graph_dim())[0]
        return self.attachment(top_level_node, dst)

    def orient(self, o):
        """ Orientation node is always labelled with -1 """
        oriented_point = copy.deepcopy(self)
        top_level_node = oriented_point.d_entities(
            oriented_point.dimension)[0]
        oriented_point.G.add_node(-1, point_class=None)
        oriented_point.G.add_edge(-1, top_level_node,
                                  edge_class=Edge(None, o=o))
        oriented_point.oriented = o
        return oriented_point
    
    def __repr__(self):
        entity_name = ["v", "e", "f", "c"]
        return entity_name[self.dimension] + str(self.id)

    def copy(self):
        return copy.deepcopy(self)



class Edge():

    def __init__(self, point, attachment=lambda x: x, o=lambda x: x):
        self.attachment = attachment
        self.point = point
        self.o = o

    def __call__(self, *x):
        return self.attachment(*self.o(x))

    def lower_dim(self):
        return self.point.dim()
    
    def __repr__(self):
        return str(self.point)
