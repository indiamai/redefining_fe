import matplotlib.pyplot as plt
import numpy as np

# temporary group implementation


def e(x): return x
# reflection on reference interval
def r(x): return -x

# reflection on triangle


def r_tri(x, y): return [-x, y]

# rotation by 60 degrees on triangle


def rot(x, y): return [-x / 2 - np.sqrt(3) * y / 2, np.sqrt(3) * x / 2 - y / 2]


class Point():

    def __init__(self, d, edges=[]):
        self.dimension = d
        if d == 0:
            assert (edges == [])
        for edge in edges:
            assert edge.lower_dim() < self.dimension
        self.edges = edges

    def dim(self):
        return self.dimension

    def hasse_diagram(self, show=True, counter=0):
        # this is currently incorrect
        if show:
            ax = plt.axes()
        plt.plot(counter, self.dimension, marker="o")
        lower_counter = counter
        for edge in self.edges:
            ax = edge.hasse_diagram(
                False, self.dimension, counter, lower_counter)
            lower_counter += 1
        if show:
            plt.show()

    def makeArrow(self, ax, mid, edge, direction=1):
        delta = 0.0001 if direction >= 0 else -0.0001
        x, y = edge(mid)
        dir_x, dir_y = edge(mid + delta)
        print("arrow")
        print((x, y))
        print((dir_x, dir_y))
        ax.arrow(x, y, x - dir_x, y - dir_y, head_width=0.05, head_length=0.1)

    def plot(self, show=True, attach=lambda x: x):
        # for now into 2 dimensional spaceo

        if self.dimension == 2:
            xs = np.linspace(-1, 1, 20)
            ax = plt.axes()
            for edge in self.edges:
                edgevals = np.array([attach(edge(x)) for x in xs])
                plt.plot(edgevals[:, 0], edgevals[:, 1])
                self.makeArrow(ax, 0, edge)
                edge.point.plot(False, edge)
        elif self.dimension == 1:
            xs = [0]
            for edge in self.edges:
                edgevals = np.array([attach(edge(x)) for x in xs])
                plt.plot(edgevals[:, 0], edgevals[:, 1], marker="o")
        else:
            raise "Error not implemented general plotting"
        plt.show()


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

    # def plot(self):


if __name__ == "__main__":
    a1 = Point(0)
    a2 = Point(0)
    a3 = Point(1, [Edge(a1, lambda x: -1), Edge(a2, lambda x: 1)])
    a4 = Point(2, [Edge(a3, lambda x: [x, -np.sqrt(3) / 3], r), Edge(a3, lambda x: [(x - 1) / 2,
               np.sqrt(3) * (3 * x + 1) / 6]), Edge(a3, lambda x: [(x + 1) / 2, np.sqrt(3) * (-3 * x + 1) / 6])])
    a4.plot()
    a4.hasse_diagram()
