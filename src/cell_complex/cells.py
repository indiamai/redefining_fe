import matplotlib.pyplot as plt
import numpy as np
import itertools
import networkx as nx


def e(x):
    # identity
    return x


def r(x):
    # reflection on reference interval
    return -x


def r_tri(xs):
    # reflection on triangle
    x, y = xs[0], xs[1]
    return [-x, y]


def rot(xs):
    # rotation by 60 degrees on triangle
    x, y = xs
    return [-x / 2 - np.sqrt(3) * y / 2, np.sqrt(3) * x / 2 - y / 2]


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
        self.G = nx.compose_all([self.G]
                                + [edge.point.get_G() for edge in edges])
        self.edges = edges

    def dim(self):
        return self.dimension

    def graph_dim(self):
        if -1 in self.G.nodes:
            dim = self.dimension + 1
        else:
            dim = self.dimension
        return dim

    def get_G(self):
        if -1 in self.G.nodes():
            temp_G = self.G.copy()
            temp_G.remove_node(-1)
            return temp_G
        return self.G

    def hasse_diagram(self, counter=0):
        ax = plt.axes()
        nx.draw_networkx(self.get_G(), pos=topo_pos(self.get_G()),
                         with_labels=True, ax=ax)
        plt.show()

    def get_d_entities(self, d, get_class=False):
        levels = [sorted(generation)
                  for generation in nx.topological_generations(self.G)]
        if get_class:
            return [self.G.nodes.data("point_class")[i]
                    for i in levels[self.graph_dim() - d]]
        return levels[self.graph_dim() - d]

    def get_dim_of_node(self, node):

        levels = [sorted(generation)
                  for generation in nx.topological_generations(self.G)]
        for i in range(len(levels)):
            if node in levels[i]:
                return self.graph_dim() - i
        raise "Error: Node not found in graph"

    def vertices(self, get_class=False):
        return self.get_d_entities(0, get_class)

    def basis_vectors(self, return_coords=False):
        vertices = self.vertices()
        top_level_node = self.get_d_entities(self.graph_dim())[0]
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
        ax.arrow(x, y, dir_x-x, dir_y-y, head_width=0.05, head_length=0.1)

    def plot(self, show=True, attach=lambda x: x):
        """ for now into 2 dimensional space """

        top_level_node = self.get_d_entities(self.graph_dim())[0]
        xs = np.linspace(-1, 1, 20)
        ax = plt.gca()
        for i in range(self.dimension - 1, -1, -1):
            nodes = self.get_d_entities(i)
            for node in nodes:
                attach = self.get_attachment_route(top_level_node, node)
                if i == 0:
                    plt.plot(attach(0)[0], attach(0)[1], 'bo')
                    plt.annotate(node, (attach(0)[0], attach(0)[1]))
                elif i == 1:
                    edgevals = np.array([attach(x) for x in xs])
                    plt.plot(edgevals[:, 0], edgevals[:, 1])
                    self.make_arrow(ax, 0, attach)
                else:
                    raise "Error not implemented general plotting"
        plt.show()

    def get_attachment_route(self, source, dst):
        paths = nx.all_simple_edge_paths(self.G, source, dst)
        attachments = [[self.G[s][d]["edge_class"]
                        for (s, d) in path] for path in paths]

        # check all attachments resolve to the same function
        source_dim = self.get_dim_of_node(source)
        basis = np.eye(source_dim)
        for i in range(source_dim):
            vals = [fold_reduce(attachment, basis[i])
                    for attachment in attachments]
            assert all(val == vals[0] for val in vals)

        return lambda x: fold_reduce(attachments[0], x)

    def orient(self, o):
        top_level_node = self.get_d_entities(self.dimension)[0]
        self.G.add_node(-1, point_class=None)
        self.G.add_edge(-1, top_level_node, edge_class=Edge(None, o=o))
        # maybe duplicate the object here rather than just modifying in place


class Edge():

    def __init__(self, point, attachment=lambda x: x, o=lambda x: x):
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
                   Edge(edges[2], lambda x: [(1 - x) / 2,
                                             np.sqrt(3) * (3 * x + 1) / 6])])
    a4.plot()
    print(a4.vertices())
    print(a4.basis_vectors(return_coords=True))
    a4.orient(rot)
    print(a4.basis_vectors(return_coords=True))
    a4.hasse_diagram()
    a4.plot()
