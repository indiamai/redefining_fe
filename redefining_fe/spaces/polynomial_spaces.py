from firedrake import *


class PolynomialSpace(object):
    """
    subdegree: the degree of the maximum degree Lagrange space that is spanned by this element. If this
    element's polynomial space does not include the constant function, this function should
    return -1.
    
    super degree: the degree of the minimum degree Lagrange space that spans this element.If this
    element contains basis functions that are not in any Lagrange space, this property should
    be None.

    Note that on a simplex cells, the polynomial space of Lagrange space is a complete polynomial
    space, but on other cells this is not true. For example, on quadrilateral cells, the degree 1
    Lagrange space includes the degree 2 polynomial xy.
    """

    def __init__(self, subdegree, superdegree):
        self.subdegree = subdegree
        self.superdegree = superdegree


    def complete(self):
        return self.subdegree == self.superdegree
    
    # def __mul__(self, x):
    #     return VectorPolynomialSpace([self, x])


class VectorPolynomialSpace(PolynomialSpace):

    def __init__(self, *spaces):
        self.component_spaces = []
        for space in spaces:
            assert isinstance(space, PolynomialSpace)
            self.component_spaces.append(space)
            # if isinstance(space, VectorPolynomialSpace):
            #     self.component_spaces.extend(space.component_spaces)

    def dim(self):
        return len(self.component_spaces)

    def complete(self):
        return all([c.complete for c in self.component_spaces])


P0 = PolynomialSpace(0, 0)
P1 = PolynomialSpace(1, 1)
P2 = PolynomialSpace(2, 2)
P3 = PolynomialSpace(3, 3)
P4 = PolynomialSpace(4, 4)

Q1 = PolynomialSpace(1, 2)
Q2 = PolynomialSpace(2, 3)
Q3 = PolynomialSpace(3, 4)
Q4 = PolynomialSpace(4, 5)