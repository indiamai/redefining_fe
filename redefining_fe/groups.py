import redefining_fe.cells as cells
from sympy.combinatorics import PermutationGroup, Permutation
from sympy.combinatorics.named_groups import SymmetricGroup, DihedralGroup, CyclicGroup, AlternatingGroup
import numpy as np
import sympy as sp


def fold_reduce(func_list, x):
    """ duplicated from cells.py """
    """ nested function composition helper function, right to left """
    prev = x
    for func in reversed(func_list):
        prev = func(prev)
    return prev


def construct_rep_func(M):
    def rep(*x):
        if isinstance(x, sp.Expr):
            breakpoint()
            x_ones = sp.r_[sp.array(*x), sp.ones(1)]
            sum = sp.matmul(x_ones, M)
        else:
            x_ones = np.r_[np.array(*x), np.ones(1)]
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
        if isinstance(x, cells.Point):
            return x.orient(self)
        # breakpoint()
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
            A = np.c_[np.array(vertices, dtype=float), np.ones(len(vertices))]
            b = np.array(vertices, dtype=float)

            M, _, _, _ = np.linalg.lstsq(A, b)
            self.identity = GroupMemberRep(base_group.identity, construct_rep_func(M), self)

            counter = 0
            for g in self.base_group.generators:
                if len(vertices) > g.size:
                    temp_perm = Permutation(g, size=len(vertices))
                    reordered = temp_perm(vertices)
                else:
                    reordered = g(vertices)
                A = np.c_[np.array(vertices, dtype=float), np.ones(len(vertices))]
                b = np.array(reordered, dtype=float)

                M, _, _, _ = np.linalg.lstsq(A, b)
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


def e(*x):
    print(isinstance(x, tuple))
    print(x)
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
