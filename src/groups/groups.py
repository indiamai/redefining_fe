import numpy as np


class Group():

    def __init__(self):
        self.members = []



def e(x):
    # identity
    return x


def r(x):
    # reflects the first component of the input array
    x[0] = -x[0]
    return x


def rot(xs, rad=2*np.pi/3):
    # rotation by rad radians
    x, y = xs[0], xs[1]
    res = [x*np.cos(rad) - y*np.sin(rad), x*np.sin(rad) + y*np.cos(rad)]
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
    def __init__(self, shape="tri"):
        super(S3, self).__init__()
        if shape == "tri":
            self.rot = rot
        elif shape == "square":
            self.rot = lambda x: rot(x, np.pi / 2)
        self.members = [e, r, self.rot,
                        lambda x: self.rot(r(x)),
                        lambda x: self.rot(self.rot(x)),
                        lambda x: self.rot(self.rot(r(x)))]


if __name__ == "__main__":
    print(rot([4, 5]))
