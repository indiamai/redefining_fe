import matplotlib.pyplot as plt
import numpy as np
import itertools
import networkx as nx

#temporary group implementation
e = lambda x: x
#reflection on reference interval
r = lambda x: -x
# reflection on triangle
r_tri = lambda x,y : [-x, y]
# rotation by 60 degrees on triangle
rot = lambda x,y : [-x/2 - np.sqrt(3)*y/2, np.sqrt(3)*x/2 - y/2]   

# helper function for hasse diagram visualisation
def topo_pos(G):
    """Display in topological order, with simple offsetting for legibility"""
    pos_dict = {}
    for i, node_list in enumerate(nx.topological_generations(G)):
        x_offset = len(node_list) / 2
        y_offset = 0
        for j, name in enumerate(node_list):
            pos_dict[name] = (j - x_offset, -i + j * y_offset)

    return pos_dict

class Point():
     
     id_iter = itertools.count()

     def __init__(self, d, edges=[]):
        self.id = next(self.id_iter)
        self.dimension = d
        if d == 0:
                assert(edges == [])

        self.G = nx.DiGraph()
        self.G.add_node(self.id, point_class = self)
        for edge in edges:
                assert edge.lower_dim() < self.dimension
                self.G.add_edge(self.id, edge.point.id, edge_class=edge)
        self.G = nx.compose_all([self.G] + [edge.point.G for edge in edges])
        self.edges = edges


     def dim(self):
        return self.dimension

     def hasse_diagram(self, show=True, counter = 0):
        if show:
                ax = plt.axes()
        nx.draw_networkx(self.G, pos = topo_pos(self.G), with_labels=True)
        if show:
                plt.show()

     def makeArrow(self, ax,mid,edge,direction=1):
             delta = 0.0001 if direction >= 0 else -0.0001
             x, y = edge(mid)
             dir_x, dir_y = edge(mid+delta)
             print("arrow")
             print((x,y))
             print((dir_x,dir_y))
             ax.arrow(x,y,x-dir_x,y-dir_y,head_width=0.05,head_length=0.1)

     def plot(self, show=True, attach= lambda x: x):
        # for now into 2 dimensional space


        if self.dimension == 2:
                xs = np.linspace(-1,1, 20)
                ax = plt.axes()
                for edge in self.edges:
                        edgevals = np.array([attach(edge(x)) for x in xs])
                        plt.plot(edgevals[:, 0], edgevals[:,1])
                        self.makeArrow(ax, 0, edge)
                        edge.point.plot(False,edge)
        elif self.dimension == 1:
                xs = [0]
                for edge in self.edges:
                        edgevals = np.array([attach(edge(x)) for x in xs])
                        plt.plot(edgevals[:, 0], edgevals[:,1], marker = "o")
        else:
                raise "Error not implemented general plotting"
        plt.show()

        
        



class Edge():

        def __init__(self, point, attachment, o = lambda x: x):
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
        edges.append(Point(1, [Edge(vertices[0], lambda x: -1), Edge(vertices[1], lambda x: 1)]))
        edges.append(Point(1, [Edge(vertices[0], lambda x: -1), Edge(vertices[2], lambda x: 1)]))
        edges.append(Point(1, [Edge(vertices[1], lambda x: -1), Edge(vertices[2], lambda x: 1)]))
        a4 = Point(2, [Edge(edges[0], lambda x: [x,-np.sqrt(3)/3]), Edge(edges[1], lambda x: [(x-1)/2, np.sqrt(3)*(3*x +1)/6]), Edge(edges[2], lambda x: [(x+1)/2, np.sqrt(3)*(-3*x +1)/6])])
        # a4.plot()
        a4.hasse_diagram()

