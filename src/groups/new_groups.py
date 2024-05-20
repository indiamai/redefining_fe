from sympy.combinatorics import PermutationGroup, Permutation
from sympy.combinatorics.named_groups import SymmetricGroup, DihedralGroup, CyclicGroup, AlternatingGroup
import numpy as np
import cell_complex.cells
import matplotlib.pyplot as plt


def fold_reduce(func_list, x):
    """ duplicated from cells.py """
    """ nested function composition helper function, right to left """
    prev = x
    for func in reversed(func_list):
        prev = func(prev)
    return prev

def construct_rep_func(M):
    def rep(*x):
        x_ones = np.r_[np.array(*x), np.ones(1)]
        # breakpoint()
        sum = np.matmul(x_ones, M)
        if len(sum.shape) > 1:
            return tuple(map(tuple, sum))
        return tuple(sum)
    return rep

class GroupMemberRep(object):

    def __init__(self, perm, rep, group):
        if not isinstance(rep, list):
            rep = [rep]

        self.perm = perm
        self.rep = rep
        self.group = group

    def __call__(self, x):
        if isinstance(x, cell_complex.cells.Point):
            return x.orient(self)
        return fold_reduce(self.rep, x)

    def permute(self, lst):
        n = len(lst)
        if n > self.perm.size:
            temp_perm = Permutation(self.perm, size=n)
            return temp_perm(lst)
        return self.perm(lst)

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

    def __init__(self, base_group, cell=None):
        assert isinstance(base_group, PermutationGroup)
        self.base_group = base_group
        self.identity = GroupMemberRep(base_group.identity, e, self)
        self.generators = []
        if cell is not None:
            self.cell = cell
            vertices = cell.vertices(return_coords=True)
            counter = 0
            for g in self.base_group.generators:
                if len(vertices) > g.size:
                    temp_perm = Permutation(g, size=len(vertices))
                    reordered = temp_perm(vertices)
                else:
                    reordered = g(vertices)
                M = np.linalg.solve(np.c_[np.array(vertices),
                                            np.ones(len(vertices))], np.array(reordered))
                rep = construct_rep_func(M)
                rep.__name__ = "g" + str(counter)
                self.generators.append(GroupMemberRep(g, rep, self))
                counter += 1

            # this order produces simpler generator lists
            self.generators.reverse()
            self._members = [self.identity]

            temp_group_elems = self.base_group._elements
            temp_group_elems.remove(self.base_group.identity)
            remaining_members = self.compute_reps(self.base_group.identity,
                                                None, temp_group_elems)
            assert (len(remaining_members) == 0)
        else:
            self.cell = None
    
    def add_cell(self, cell):
        return GroupRepresentation(self.base_group, cell=cell)
        
    def members(self, perm=False):
        if self.cell is None:
            raise ValueError("Group does not have a domain - members have not been calculated")
        if perm:
            return [m.perm for m in self._members]
        return self._members

    def size(self):
        assert len(self._members) == self.base_group.order()
        return self.base_group.order()
    
    def transform_between_perms(self, perm1, perm2):
        member_perms = self.members(perm=True)
        # breakpoint()
        perm1 = Permutation.from_sequence(perm1)
        perm2 = Permutation.from_sequence(perm2)
        assert perm1 in member_perms
        assert perm2 in member_perms
        # assert list(perm2) in [m.array_form for m in member_perms]
        return self.get_member(~Permutation(perm1)) * self.get_member(Permutation(perm2))

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
                    self._members.append(generator)
                else:
                    new_path = path.copy()
                    new_path.extend(generator.rep)
                    self._members.append(GroupMemberRep(new_perm,
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
        for m in self.members():
            if m.perm == perm:
                return m
        raise ValueError("Permutation not a member of group")
    
    def __mul__(self, other_group):
        return GroupRepresentation(PermutationGroup(self.base_group.generators + other_group.base_group.generators))

    def __truediv__(self, other_frac):
        """ This isn't a mathematically accurate representation of
            what it means to be a quotient group but it works on S3/S2
            Have to compare cyclic forms as groups may not be defined on
            the same number of elements
            Doesn't work on D4/S2 but does on D4/C4 """
        assert isinstance(other_frac, GroupRepresentation)
        self_cyclic_gens = [gen.cyclic_form
                            for gen in self.base_group.generators]
        other_cyclic_gens = [gen.cyclic_form
                             for gen in other_frac.base_group.generators]
        # # breakpoint()
        # print(all([c2 in self_cyclic_gens for c2 in other_cyclic_gens]))
        if not all([c2 in self_cyclic_gens for c2 in other_cyclic_gens]):
            raise ValueError("Invalid Quotient - mismatched cycles")
        remaining_perms = [gen for gen in self.base_group.generators
                           if gen.cyclic_form not in other_cyclic_gens]

        if len(remaining_perms) == 0:
            raise ValueError("Invalid Quotient - no group formed")

        return GroupRepresentation(PermutationGroup(remaining_perms))

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

def r_y(x):
    # reflection in the second component
    if isinstance(x, int) or len(x) < 2:
        raise ValueError("Input array too short")
    x_list = list(x)
    x_list[1] = -x_list[1]
    return tuple(x_list)

def rot(xs, rad=2*np.pi/3):
    #  anticlockwise rotation by rad radians, default is 120 deg
    x, y = xs[0], xs[1]
    res = (x*np.cos(rad) - y*np.sin(rad), x*np.sin(rad) + y*np.cos(rad))
    return res


def sqrot(xs):
    # 90 degree rotation
    return rot(xs, np.pi / 2)


def g1(xs):
    raise NotImplementedError("Tetrahedron implementation incomplete")


def g2(xs):
    # 120 degree rotation clockwise
    raise NotImplementedError("Tetrahedron implementation incomplete")
    return rot(xs, - 2*np.pi / 3)


S1 = GroupRepresentation(SymmetricGroup(1))
S2 = GroupRepresentation(SymmetricGroup(2))
S3 = GroupRepresentation(SymmetricGroup(3))
S4 = GroupRepresentation(SymmetricGroup(4))

D4 = GroupRepresentation(DihedralGroup(4))

C3 = GroupRepresentation(CyclicGroup(3))
C4 = GroupRepresentation(CyclicGroup(4))

Z2 = GroupRepresentation(CyclicGroup(2))
Z3 = GroupRepresentation(CyclicGroup(3))
Z4 = GroupRepresentation(CyclicGroup(4))


D2 = GroupRepresentation(DihedralGroup(2))
A4 = GroupRepresentation(AlternatingGroup(4))
A3 = GroupRepresentation(AlternatingGroup(3))

 

# S2 = GroupRepresentation(SymmetricGroup(2), {Permutation(0, 1): r})
# S3 = GroupRepresentation(SymmetricGroup(3), {Permutation(0, 1, 2): rot,
#                                              Permutation([1, 0, 2]): r})

# D4 = GroupRepresentation(DihedralGroup(4), {Permutation(0, 1, 2, 3): sqrot,
#                                             Permutation([3, 2, 1, 0]): r})

# C3 = GroupRepresentation(CyclicGroup(3), {Permutation(0, 1, 2): rot})
# C4 = GroupRepresentation(CyclicGroup(4), {Permutation(0, 1, 2, 3): sqrot})

# S4 = GroupRepresentation(SymmetricGroup(4), {Permutation(0, 1, 2, 3): g1,
#                                              Permutation([1, 0, 2, 3]): g2})

# if __name__ == "__main__":

    
    # print(C3.members)
    # print(C4.members)
    # print((D4/C4).members)
    # print((S3/S2).base_group._elements)

    # for m in S3.members:
    #     print(m)
    #     print(m((-1, -np.sqrt(3)/3)))
    #     coord = m((-0.9, -np.sqrt(3)/3))
    #     plt.scatter(coord[0], coord[1], marker="o")
    #     plt.text(coord[0], coord[1], repr(m))

    # plt.show()

    # source = (-1, -0.9)
    # for m in D4.members:
    #     print(m)
    #     print(m(source))
    #     coord = m(source)
    #     plt.scatter(coord[0], coord[1], marker="o")
    #     plt.text(coord[0], coord[1], repr(m))
    # plt.plot([-1, -0.5, 0, 0.5, 1], [-1, -1, -1, -1, -1])
    # plt.plot([-1, -0.5, 0, 0.5, 1], [1, 1, 1, 1, 1])
    # plt.plot([-1, -1, -1, -1, -1], [-1, -0.5, 0, 0.5, 1])
    # plt.plot([1, 1, 1, 1, 1], [-1, -0.5, 0, 0.5, 1])

    # plt.show()

    # gens = refD4.generators
    # d4rot = gens[1]
    # d4r = gens[0]

    # point = d4rot * d4rot * d4rot * d4rot
    # print(point)
    # print(point(coord))

# a = Permutation(0, 1)(2, 3)
# b = Permutation(0, 2)(1, 3)
# c = Permutation(0, 3)(1, 2)
# K4 = GroupRepresentation(PermutationGroup(a, b, c))
# print("K4")
# print(K4.base_group.generators)
# print(K4.base_group.order())

# print("S4/S2")
# print((S4/S2).base_group.generators)
# print((S4/S2).base_group.order())
# c3_2_tet = (S4/S2).add_cell(tetrahedron)
# mems = c3_2_tet.members()
# for m in mems:
#     print(m.perm.array_form)
# # # breakpoint()

# a = Permutation(1)(0, 2, 3)
# b = Permutation(2, 3)
# C3_2 = GroupRepresentation(PermutationGroup(a))
# print("Weird group")
# print(C3_2.base_group.generators)
# print(C3_2.base_group.order())
# c3_2_tet = C3_2.add_cell(tetrahedron)
# mems = c3_2_tet.members()
# for m in mems:
#     print(m.perm.array_form)
# # # breakpoint()
# # print((D2 * C3).base_group.generators)

