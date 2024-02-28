from sympy.combinatorics import PermutationGroup, Permutation
from sympy.combinatorics.named_groups import SymmetricGroup, DihedralGroup, CyclicGroup
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

    def __init__(self, base_group, reps_dict):
        assert isinstance(base_group, PermutationGroup)
        self.base_group = base_group
        self.identity = GroupMemberRep(base_group.identity, e, self)
        self.generators = []

        for perm in reps_dict.keys():
            if perm in base_group.generators:
                self.generators.append(GroupMemberRep(perm,
                                                      reps_dict[perm], self))
            else:
                raise ValueError("Generator key does not match any generator in the base group")

        # this order produces simpler generator lists
        self.generators.reverse()
        self.members = [self.identity]

        temp_group_elems = base_group._elements
        temp_group_elems.remove(base_group.identity)
        remaining_members = self.compute_reps(base_group.identity,
                                              None, temp_group_elems)
        assert (len(remaining_members) == 0)

    def size(self):
        assert len(self.members) == self.base_group.order()
        return self.base_group.order()

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

    def __truediv__(self, other_frac):
        """ This isn't a mathematically accurate representation of
            what it means to be a quotient group but it works on S3/S2
            Have to compare cyclic forms as groups may not be defined on
            the same number of elements
            Doesn't work on D4/S2 but does on D4/C4 """
        assert isinstance(other_frac, GroupRepresentation)
        self_cyclic_gens = [gen.perm.cyclic_form
                            for gen in self.generators]
        other_cyclic_gens = [gen.perm.cyclic_form
                             for gen in other_frac.generators]

        if not all([c2 in self_cyclic_gens for c2 in other_cyclic_gens]):
            raise ValueError("Invalid Quotient - mismatched cycles")
        remaining_perms = [gen.perm for gen in self.generators
                           if gen.perm.cyclic_form not in other_cyclic_gens]
        rep_dict = {gen.perm: gen.rep for gen in self.generators
                    if gen.perm.cyclic_form not in other_cyclic_gens}

        if len(remaining_perms) == 0:
            raise ValueError("Invalid Quotient - no group formed")

        return GroupRepresentation(PermutationGroup(remaining_perms),
                                   rep_dict)

# Function Representation of the coordinate transforms that make up the groups.


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
    # rotation by rad radians, default is 60 deg
    x, y = xs[0], xs[1]
    res = (x*np.cos(rad) - y*np.sin(rad), x*np.sin(rad) + y*np.cos(rad))
    return res


def sqrot(xs):
    # 90 degree rotation
    return rot(xs, np.pi / 2)


S1 = GroupRepresentation(SymmetricGroup(1), {})
S2 = GroupRepresentation(SymmetricGroup(2), {Permutation(0, 1): r})
S3 = GroupRepresentation(SymmetricGroup(3), {Permutation(0, 1, 2): rot,
                                             Permutation([1, 0, 2]): r})

D4 = GroupRepresentation(DihedralGroup(4), {Permutation(0, 1, 2, 3): sqrot,
                                            Permutation([3, 2, 1, 0]): r})

C3 = GroupRepresentation(CyclicGroup(3), {Permutation(0, 1, 2): rot})
C4 = GroupRepresentation(CyclicGroup(4), {Permutation(0, 1, 2, 3): sqrot})

if __name__ == "__main__":
    # print(C3.members)
    # print(C4.members)
    # print((D4/C4).members)
    # print((S3/S2).base_group._elements)

    for m in S3.members:
        print(m)
        print(m((-1, -np.sqrt(3)/3)))
        coord = m((-0.9, -np.sqrt(3)/3))
        plt.scatter(coord[0], coord[1], marker="o")
        plt.text(coord[0], coord[1], repr(m))

    plt.show()

    source = (-1, -0.9)
    for m in D4.members:
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

    # gens = refD4.generators
    # d4rot = gens[1]
    # d4r = gens[0]

    # point = d4rot * d4rot * d4rot * d4rot
    # print(point)
    # print(point(coord))
