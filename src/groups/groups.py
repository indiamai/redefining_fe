import numpy as np


class Group():

    def __init__(self, members=[]):
        self.members = members


class GroupMember():
    def __init__(self, name):
        self.name = name


# class e(GroupMember):
#     # identity

#     def __call__(x):
#         return x

def e(x):
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


class clsS1(Group):

    def __init__(self):
        super(clsS1, self).__init__()
        self.members = [e]
        

class clsS2(Group):
    def __init__(self):
        super(clsS2, self).__init__()
        self.members = [e, r]


class clsS3(Group):

    def __init__(self):
        self.rot = rot
        super(clsS3, self).__init__(members=[e, r, self.rot,
                                             lambda x: self.rot(r(x)),
                                             lambda x: self.rot(self.rot(x)),
                                             lambda x: self.rot(self.rot(r(x)))])

    def __truediv__(self, other_frac):
        if isinstance(other_frac, clsS2):
            return Group([e, self.rot, lambda x: self.rot(self.rot(x))])
        else:
            raise "Unknown quotient"


class clsD4(Group):
    def __init__(self):
        self.rot = lambda x: rot(x, np.pi / 2)
        super(clsD4, self).__init__(members=[e, r, self.rot,
                                             lambda x: self.rot(r(x)),
                                             lambda x: self.rot(self.rot(x)),
                                             lambda x: self.rot(self.rot(r(x))),
                                             lambda x: self.rot(
                                                       self.rot(self.rot(x))),
                                             lambda x: self.rot(
                                                       self.rot(self.rot(r(x))))])

    def __truediv__(self, other_frac):
        if isinstance(other_frac, clsS2):
            return Group([e, self.rot,
                          lambda x: self.rot(self.rot(x)),
                          lambda x: self.rot(self.rot(self.rot(x)))])
        else:
            raise "Unknown quotient"


S1 = clsS1()
S2 = clsS2()
S3 = clsS3()
D4 = clsD4()
if __name__ == "__main__":
    print(rot([4, 5]))
