import redefining_fe.cells as cells
from sympy.combinatorics import PermutationGroup, Permutation
from sympy.combinatorics.named_groups import SymmetricGroup, DihedralGroup, CyclicGroup, AlternatingGroup
import numpy as np
import sympy as sp
import math
from redefining_fe.utils import fold_reduce


def construct_rep_func(M):
    """
    Convert a matrix to a representation function

    Args:
        M: a numpy matrix

    Returns a function that applies the matrix to the arguments (which can be symbolic or numeric) """
    def rep(*x):
        if isinstance(x, sp.Expr):
            x_ones = sp.r_[sp.array(x), sp.ones(1)]
            sum = sp.matmul(x_ones, M)
        else:
            x_ones = np.r_[np.array(x), np.ones(1)]
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
        return fold_reduce(self.rep, *x)

    def permute(self, lst):
        n = len(lst)
        if n > self.perm.size:
            temp_perm = Permutation(self.perm, size=n)
            return temp_perm(lst)
        return self.perm(lst)

    def compute_num_rep(self, val_list=None):
        m_array = self.perm.array_form
        identity = self.group.identity.perm.array_form
        if not val_list:
            val_list = self.perm.array_form
        else:
            val_list = self.permute(val_list)
        val = 0
        for i in range(len(identity)):
            loc = m_array.index(identity[i])
            m_array.remove(identity[i])
            val += loc * math.factorial(len(identity) - i - 1)
        return val, val_list

    def __mul__(self, x):
        assert isinstance(x, GroupMemberRep)
        return self.group.get_member(self.perm * x.perm)

    def __repr__(self):
        string = ""
        for rep in self.rep:
            string += rep.__name__
            string += " "
        return str(self.perm.array_form)

    def __eq__(self, value):
        assert isinstance(value, GroupMemberRep)
        return self.perm == value.perm

    # def _to_dict(self):
    #     o_dict = {self.dict_id(): self.perm}
    #     return o_dict

    # def dict_id(self):
    #     return "GroupMemberRep" + str(id(self))


class GroupRepresentation(object):
    """
    A representation of a group by its matrix operations.

    Args:
        base_group: the sympy group that is being represented
        cell (optional): the cell the group is representing the operations on

    """

    def __init__(self, base_group, cell=None):
        assert isinstance(base_group, PermutationGroup)
        self.base_group = base_group
        self.identity = GroupMemberRep(base_group.identity, e, self)
        self.generators = []
        if cell is not None:
            self.cell = cell
            vertices = cell.vertices(return_coords=True)
            try:
                A = np.c_[np.array(vertices, dtype=float), np.ones(len(vertices))]
            except ValueError:
                breakpoint()
            b = np.array(vertices, dtype=float)

            M, _, _, _ = np.linalg.lstsq(A, b, rcond=None)
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

                M, _, _, _ = np.linalg.lstsq(A, b, rcond=None)
                rep = construct_rep_func(M)
                rep.__name__ = "g" + str(counter)
                self.generators.append(GroupMemberRep(g, rep, self))
                counter += 1

            # this order produces simpler generator lists
            self.generators.reverse()
            self._members = [self.identity]

            temp_group_elems = self.base_group.elements

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

    def compute_num_reps(self, base_val=0):
        """ Computed the numerical represention of each member as compared to the identity.
        Where the numerical rep is:

        M.index(id[0]) = a; M.remove(id[0])
        M.index(id[1]) = b; M.remove(id[1])
        M.index(id[2]) = c; M.remove(id[2])

        o = (a * 2!) + (b * 1!) + (c * 0!)
        """
        members = self.members()
        res = {}
        for m in members:
            val, perm = m.compute_num_rep(base_val)
            res[val] = perm
        return res

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
        if not all([c2 in self_cyclic_gens for c2 in other_cyclic_gens]):
            raise ValueError("Invalid Quotient - mismatched cycles")
        remaining_perms = [gen for gen in self.base_group.generators
                           if gen.cyclic_form not in other_cyclic_gens]

        if len(remaining_perms) == 0:
            raise ValueError("Invalid Quotient - no group formed")

        return GroupRepresentation(PermutationGroup(remaining_perms))

    def __repr__(self):
        return "GR"

    # def __eq__(self, other):
    #     # TODO work on idea of group equality
    #     assert isinstance(other, GroupRepresentation)
    #     res = True
    #     for m in self.members():
    #         res = res and m in other.members()
    #     return res

    def _to_dict(self):
        return {"members": [m.perm.array_form for m in self._members]}

    def dict_id(self):
        return "Group"

    def _from_dict(o_dict):
        perm_group = PermutationGroup([Permutation(m) for m in o_dict["members"]])
        # , o_dict["cell"]
        return GroupRepresentation(perm_group)


# Function Representation of the coordinate transforms that make up the groups.


def e(*x):
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


def get_sym_group(n):
    return GroupRepresentation(SymmetricGroup(n))


def get_cyc_group(n):
    return GroupRepresentation(CyclicGroup(n))


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
