import matplotlib.pyplot as plt
import numpy as np
import itertools
import networkx as nx
from groups.groups import e, r, rot
import copy


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
        prev = func(*prev)
    return prev


def normalise(v):
    norm = np.linalg.norm(v)
    return v / norm


def make_arrow(ax, mid, edge, direction=1):
    delta = 0.0001 if direction >= 0 else -0.0001
    x, y = edge(mid)
    dir_x, dir_y = edge(mid + delta)
    ax.arrow(x, y, dir_x-x, dir_y-y, head_width=0.05, head_length=0.1)


class Point():

    id_iter = itertools.count()

    def __init__(self, d, edges=[], oriented=False):
        self.id = next(self.id_iter)
        self.dimension = d
        if d == 0:
            assert (edges == [])

        self.oriented = oriented

        self.G = nx.DiGraph()
        self.G.add_node(self.id, point_class=self)
        for edge in edges:
            assert edge.lower_dim() < self.dimension
            self.G.add_edge(self.id, edge.point.id, edge_class=edge)
        self.G = nx.compose_all([self.G]
                                + [edge.point.graph() for edge in edges])
        self.connections = edges

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
    
    def vertices(self, get_class=False):
        if self.oriented:
            return self.oriented.permute(self.d_entities(0, get_class))
        return self.d_entities(0, get_class)
    
    def edges(self, get_class=False):
        if self.oriented:
            return self.oriented.permute(self.d_entities(1, get_class))
        return self.d_entities(1, get_class)

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


    def plot(self, show=True, plain=False):
        """ for now into 2 dimensional space """

        top_level_node = self.d_entities(self.graph_dim())[0]
        xs = np.linspace(-1, 1, 20)
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
            for node in nodes:
                attach = self.attachment(top_level_node, node)
                if i == 0:
                    plotted = attach()
                    if len(plotted) < 2:
                        plotted = (plotted[0], 0)
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
                elif i == 2:
                    # write surface plot here
                    continue
                else:
                    raise "Error not implemented general plotting"
        if show:
            plt.show()

    def plot3d(self):
        # only plots vertices
        assert self.dimension == 3
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')

        top_level_node = self.d_entities(self.graph_dim())[0]
        nodes = self.d_entities(0)
        for node in nodes:
            attach = self.attachment(top_level_node, node)
            plotted = attach()
            print(plotted)
            ax.scatter(plotted[0], plotted[1], plotted[2])
        plt.show()

    def attachment(self, source, dst):
        paths = nx.all_simple_edge_paths(self.G, source, dst)
        attachments = [[self.G[s][d]["edge_class"]
                        for (s, d) in path] for path in paths]

        if len(attachments) == 0:
            raise ValueError("No paths from node {} to node {}"
                             .format(source, dst))
        
        # check all attachments resolve to the same function
        if len(attachments) > 1:
            dst_dim = self.dim_of_node(dst)
            basis = np.eye(dst_dim)
            if dst_dim == 0:
                vals = [fold_reduce(attachment) for attachment in attachments]
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

    def __init__(self, point, attachment=e, o=e):
        self.attachment = attachment
        self.point = point
        self.o = o

    def __call__(self, *x):
        return self.attachment(*self.o(x))

    def lower_dim(self):
        return self.point.dim()
    
    def __repr__(self):
        return str(self.point)
