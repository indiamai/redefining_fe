import matplotlib.pyplot as plt
import numpy as np

class Point():

     def __init__(self, d, edges=[]):
        self.dimension = d
        if d == 0:
                assert(edges == [])
        for edge in edges:
                assert edge.lower_dim() < self.dimension
        self.edges = edges


     def dim(self):
        return self.dimension

     def hasse_diagram(self, show=True, counter = 0):
        if show:
                ax = plt.figure()
        plt.plot(self.dimension, counter, marker = "o")
        lower_counter = 0
        for edge in self.edges:
                ax = edge.hasse_diagram(False, self.dimension, counter, lower_counter)
                lower_counter += 1
        if show:
                plt.show()

     def plot(self):
        # for now into 2 dimensional space


        if self.dimension == 2:
                xs = np.linspace(-1,1, 20)
                for edge in self.edges:
                        edgevals = np.array([edge.attachment(x) for x in xs])
                        plt.plot(edgevals[:, 0], edgevals[:,1])
        else:
                raise "Error not implemented general plotting"
        plt.show()

        
        



class Edge():

        def __init__(self, point, attachment, o = "identity"):
                self.attachment = attachment
                self.point = point
                self.o = o

        def lower_dim(self):
                return self.point.dim()

        def hasse_diagram(self, show, dim, upper_counter, lower_counter):
                plt.plot([upper_counter, lower_counter], [dim, self.point.dim()])
                return self.point.hasse_diagram(show, lower_counter)

        def plot(self

if __name__ == "__main__":
        a1 = Point(0)
        a2 = Point(0)
        a3 = Point(1, [Edge(a1, lambda x: -1), Edge(a2, lambda x: 1)])
        a4 = Point(2, [Edge(a3, lambda x: [x,-np.sqrt(3)/3]), Edge(a3, lambda x: [(x-1)/2, np.sqrt(3)*(3*x +1)/6]), Edge(a3, lambda x: [(x+1)/2, np.sqrt(3)*(-3*x +1)/6])])
        a4.plot()
#        a3.hasse_diagram()


