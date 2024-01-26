import matplotlib.pyplot as plt
import numpy as np
import itertools
import networkx as nx

# temporary group implementation


def e(x): return x

# reflection on reference interval


def r(x): return -x

# reflection on triangle


def r_tri(x, y): return [-x, y]

# rotation by 60 degrees on triangle


def rot(x, y): return [-x / 2 - np.sqrt(3) * y / 2, np.sqrt(3) * x / 2 - y / 2]

# helper function for hasse diagram visualisation


def topo_pos(G):
    """Display in topological order, with simple offsetting for legibility"""
    pos_dict = {}
    for i, node_list in enumerate(nx.topological_generations(G)):
        x_offset = len(node_list) / 2
        y_offset = 0
        for j, name in enumerate(node_list):
            pos_dict[name] = (j - x_offset, i - j * y_offset)

    return pos_dict


def fold_reduce(func_list, x):
    """ nested function composition helper function, right to left """
    prev = x
    for func in reversed(func_list):
        prev = func(prev)
    return prev


class Point():

    id_iter = itertools.count()

    def __init__(self, d, edges=[]):
        self.id = next(self.id_iter)
        self.dimension = d
        if d == 0:
            assert (edges == [])

        self.G = nx.DiGraph()
        self.G.add_node(self.id, point_class=self)
        for edge in edges:
            assert edge.lower_dim() < self.dimension
            self.G.add_edge(self.id, edge.point.id, edge_class=edge)
        self.G = nx.compose_all([self.G] + [edge.point.G for edge in edges])
        self.edges = edges

    def dim(self):
        return self.dimension

    def hasse_diagram(self, counter=0):
        ax = plt.axes()
        nx.draw_networkx(self.G, pos=topo_pos(self.G), with_labels=True, ax=ax)
        plt.show()

    def get_d_entities(self, d, get_class=False):
        levels = [sorted(generation)
                  for generation in nx.topological_generations(self.G)]
        if get_class:
            return [self.G.nodes.data("point_class")[i]
                    for i in levels[self.dimension - d]]
        return levels[self.dimension - d]

    def vertices(self, get_class=False):
        return self.get_d_entities(0, get_class)

    def basis_vectors(self, return_coords=False):
        vertices = self.vertices()
        top_level_node = self.get_d_entities(self.dimension)[0]
        v_0 = vertices[0]
        if return_coords:
            v_0_coords = np.array(
                (self.get_attachment_route(top_level_node, v_0)(0)))
        basis_vecs = []
        for v in vertices[1:]:
            if return_coords:
                v_coords = np.array(
                    (self.get_attachment_route(top_level_node, v)(0)))
                basis_vecs.append(v_coords - v_0_coords)
            else:
                basis_vecs.append((v, v_0))
        return basis_vecs

    def make_arrow(self, ax, mid, edge, direction=1):
        delta = 0.0001 if direction >= 0 else -0.0001
        x, y = edge(mid)
        dir_x, dir_y = edge(mid + delta)
        ax.arrow(x, y, x - dir_x, y - dir_y, head_width=0.05, head_length=0.1)

    def plot(self, show=True, attach=lambda x: x):
        """ for now into 2 dimensional space """

        top_level_node = self.get_d_entities(self.dimension)[0]
        xs = np.linspace(-1, 1, 20)

        for i in range(self.dimension - 1, -1, -1):
            nodes = self.get_d_entities(i)
            for node in nodes:
                attach = self.get_attachment_route(top_level_node, node)
                if i == 0:
                    plt.plot(attach(0)[0], attach(0)[1], 'bo')
                elif i == 1:
                    edgevals = np.array([attach(x) for x in xs])
                    plt.plot(edgevals[:, 0], edgevals[:, 1])
                    # this is missing arrows on edges
                else:
                    raise "Error not implemented general plotting"
        plt.show()

    def get_attachment_route(self, source, dst):
        # add assertion that the paths are the same
        paths = nx.all_simple_edge_paths(self.G, source, dst)
        attachments = [[self.G[s][d]["edge_class"]
                        for (s, d) in path] for path in paths]
        return lambda x: fold_reduce(attachments[0], x)

    def orient(self, o):
        top_level_node = self.get_d_entities(self.dimension)[0]
        print(self.G.edges(top_level_node))


class Edge():

    def __init__(self, point, attachment, o=lambda x: x):
        self.attachment = attachment
        self.point = point
        self.o = o

    def __call__(self, x):
        return self.attachment(self.o(x))

    def lower_dim(self):
        return self.point.dim()

    def hasse_diagram(self, show, dim, upper_counter, lower_counter):
        plt.plot([upper_counter, lower_counter], [dim, self.point.dim()])
        return self.point.hasse_diagram(show, lower_counter)


if __name__ == "__main__":
    vertices = []
    for i in range(3):
        vertices.append(Point(0))
    edges = []
    edges.append(
        Point(1, [Edge(vertices[0], lambda x: -1),
                  Edge(vertices[1], lambda x: 1)]))
    edges.append(
        Point(1, [Edge(vertices[0], lambda x: -1),
                  Edge(vertices[2], lambda x: 1)]))
    edges.append(
        Point(1, [Edge(vertices[1], lambda x: -1),
                  Edge(vertices[2], lambda x: 1)]))

    a4 = Point(2, [Edge(edges[0], lambda x: [x, -np.sqrt(3) / 3]),
                   Edge(edges[1], lambda x: [(x - 1) / 2,
                                             np.sqrt(3) * (3 * x + 1) / 6]),
                   Edge(edges[2], lambda x: [(x + 1) / 2,
                                             np.sqrt(3) * (-3 * x + 1) / 6])])
#     a4.plot()
    # print(a4.vertices())
    # print(a4.basis_vectors())
    # print(a4.basis_vectors(return_coords=True))
    a4.hasse_diagram()
#     a4.orient(lambda x: x)
