from sympy.combinatorics import PermutationGroup
from sympy.combinatorics.named_groups import SymmetricGroup, DihedralGroup
import numpy as np
import matplotlib.pyplot as plt


def fold_reduce(func_list, x):
    """ duplicated from cells.py """
    """ nested function composition helper function, right to left """
    prev = x
    for func in reversed(func_list):
        prev = func(prev)
    return prev


class GroupMemberRep(object):

    def __init__(self, perm, rep, group):
        if not isinstance(rep, list):
            rep = [rep]

        self.perm = perm
        self.rep = rep
        self.group = group

    def __call__(self, x):
        return fold_reduce(self.rep, x)

    def __mul__(self, x):
        assert isinstance(x, GroupMemberRep)
        return self.group.get_member(self.perm * x.perm)

    def __repr__(self):
        string = ""
        for rep in self.rep:
            string += rep.__name__
            string += " "
        return string


class GroupRepresentation(object):

    def __init__(self, base_group, generator_reps):
        assert isinstance(base_group, PermutationGroup)
        self.identity = GroupMemberRep(base_group.identity, e, self)
        self.generators = []
        for (perm, rep) in zip(base_group.generators, generator_reps):
            self.generators.append(GroupMemberRep(perm, rep, self))

        # this order produces simpler generator lists
        self.generators.reverse()
        self.members = [self.identity]

        temp_group_elems = base_group._elements
        temp_group_elems.remove(base_group.identity)
        remaining_members = self.compute_reps(base_group.identity, None, temp_group_elems)
        assert (len(remaining_members) == 0)

        print(self.members)

    def compute_reps(self, g, path, remaining_members):
        # breadth first search to find generator representations of all members
        if len(remaining_members) == 0:
            return remaining_members

        next_candidates = []
        for generator in self.generators:
            new_perm = g*generator.perm
            if new_perm in remaining_members:
                if not path:
                    new_path = generator.rep
                    assert (new_perm == generator.perm)
                    self.members.append(generator)
                else:
                    new_path = path.copy()
                    new_path.extend(generator.rep)
                    self.members.append(GroupMemberRep(new_perm,
                                                       new_path,
                                                       self))
                remaining_members.remove(new_perm)
                next_candidates.append((new_perm, new_path))

        for (new_perm, new_path) in next_candidates:
            remaining_members = self.compute_reps(new_perm,
                                                  new_path,
                                                  remaining_members)
        return remaining_members

    def get_member(self, perm):
        for m in self.members:
            if m.perm == perm:
                return m
        raise ValueError("Permutation not a member of group")


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


def sqrot(xs):
    return rot(xs, np.pi / 2)


if __name__ == "__main__":
    s3 = SymmetricGroup(3)
    myS3 = GroupRepresentation(s3, [rot, r])

    for m in myS3.members:
        print(m)
        print(m((-1, -np.sqrt(3)/3)))
        coord = m((-0.9, -np.sqrt(3)/3))
        plt.scatter(coord[0], coord[1], marker="o")
        plt.text(coord[0], coord[1], repr(m))

    plt.show()

    d4 = DihedralGroup(4)
    myD4 = GroupRepresentation(d4, [sqrot, r])
    source = (-1, -0.9)
    for m in myD4.members:
        print(m)
        print(m(source))
        coord = m(source)
        plt.scatter(coord[0], coord[1], marker="o")
        plt.text(coord[0], coord[1], repr(m))
    plt.plot([-1, -0.5, 0, 0.5, 1], [-1, -1, -1, -1, -1])
    plt.plot([-1, -0.5, 0, 0.5, 1], [1, 1, 1, 1, 1])
    plt.plot([-1, -1, -1, -1, -1], [-1, -0.5, 0, 0.5, 1])
    plt.plot([1, 1, 1, 1, 1], [-1, -0.5, 0, 0.5, 1])

    plt.show()

    gens = myD4.generators
    myrot = gens[1]
    myr = gens[0]
    # rotr = myrot * myr
    # print("rotr")
    # print(rotr.perm)
    # print(rotr)

    point = myrot * myrot * myrot * myrot
    print(point.perm)
    print(point)
    print(point(coord))
