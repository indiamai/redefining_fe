from src.groups.group import *


class Element():

    def __init__(self, g, p):
        self.group = g
        self.permutation = p

    def __mul__(self, e):
        assert isinstance(e, self.__class__)
        assert e.group.equals(self.group)
        return self.__class__(self.group, self.permutation * e.permutation)

    def __repr__(self):
        return repr(self.permutation)


class Matrix_(Element):

    def __init__(self, g, p, m):
        self.matrix = m
        super(self, Matrix_Representation).__init__(g, p)


if __name__ == "__main__":
    g = S2()
    print(g.elements())
    p = Permutation([1, 0])
    r = Element(g, p)
    a = r * r
    print("a", a)
    print("perms")
    e = Permutation([0, 1])
    print(e)
    print(p)
    print(p * e)
    print(p * p)
