import numpy as np


class Group():

    def __init__(self, members=[]):
        self.members = members


def e(x):
    # identity
    return x


def r(x):
    # reflection in the first component
    if isinstance(x, int):
        x = [x]
    x_list = list(x)
    x_list[0] = -x_list[0]
    return tuple(x_list)


def rot(xs, rad=2*np.pi/3):
    # rotation by rad radians
    x, y = xs[0], xs[1]
    res = (x*np.cos(rad) - y*np.sin(rad), x*np.sin(rad) + y*np.cos(rad))
    return res


class S1(Group):

    def __init__(self):
        super(S1, self).__init__()
        self.members = [lambda x: x]


class S2(Group):
    def __init__(self):
        super(S2, self).__init__()
        self.members = [e, r]


class S3(Group):

    def __init__(self):
        self.rot = rot
        super(S3, self).__init__(members=[e, r, self.rot,
                                          lambda x: self.rot(r(x)),
                                          lambda x: self.rot(self.rot(x)),
                                          lambda x: self.rot(self.rot(r(x)))])

    def __truediv__(self, other_frac):
        if isinstance(other_frac, S2):
            return Group([e, self.rot, lambda x: self.rot(self.rot(x))])
        else:
            raise "Unknown quotient"


class D4(Group):
    def __init__(self):
        self.rot = lambda x: rot(x, np.pi / 2)
        super(D4, self).__init__(members=[e, r, self.rot,
                                          lambda x: self.rot(r(x)),
                                          lambda x: self.rot(self.rot(x)),
                                          lambda x: self.rot(self.rot(r(x))),
                                          lambda x: self.rot(
                                                    self.rot(self.rot(x))),
                                          lambda x: self.rot(
                                                    self.rot(self.rot(r(x))))])

    def __truediv__(self, other_frac):
        if isinstance(other_frac, S2):
            return Group([e, self.rot,
                          lambda x: self.rot(self.rot(x)),
                          lambda x: self.rot(self.rot(self.rot(x)))])
        else:
            raise "Unknown quotient"


if __name__ == "__main__":
    print(rot([4, 5]))
