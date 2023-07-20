from sympy.combinatorics.named_groups import SymmetricGroup
from sympy.combinatorics import Permutation, PermutationGroup

class Group():
	
	def __init__(self, G):
		assert G.is_group
		self.G = G
	
	def quotient(self, N):
		print("Quotient")
		print(g.G.elements)
		print(N.G.elements)
		print(self.G.coset_transversal(N.G))

	def elements(self):
		return self.G.elements

	def is_element(self, e):
		return self.G.contains(e)


class S1(Group):
	def __init__(self):
		super(S1, self).__init__(SymmetricGroup(1))

	def perm_name(self, p):
		assert self.G.contains(p)
		if p == self.G.identity:
			return "e"

class S2(Group):
	def __init__(self):
		super(S2, self).__init__(SymmetricGroup(2))

	def perm_name(self, p):
		assert self.G.contains(p)
		if p == self.G.identity:
			return "e"
		else:
			return "r"
class S3(Group):
	def __init__(self):
		super(S3, self).__init__(SymmetricGroup(3))

	def perm_name(self, p):
		assert self.G.contains(p)
		rank_list = {0: "e", 1: "r", 2:"rot^2r", 3: "rot^2", 4: "rot", 5: "rotr"}
		return rank_list[p.rank()]

class S4(Group):	

	def __init__(self):
		super(S4, self).__init__(SymmetricGroup(4))


if __name__== "__main__":
	g = S3() 
	print(g.G.is_group)
	print(g.perm_name(Permutation([0,2,1])))
	exit()
	print(g.G.conjugacy_classes())
	print("elements",g.G.elements)
	print("id",g.G.identity)
	print("gen",g.G.generators)
	g2 = Group(SymmetricGroup(2))
	print(g2.G.is_subgroup(g.G, strict = False))
	p = Permutation([0,2,1])
	print(g.G.contains(p))
	g3 = g.G.subgroup([p])
	print(g.G.coset_factor(p))
	
	perm = PermutationGroup([p])
	N= Group(perm)
	g.quotient(N)
	tet = Group(SymmetricGroup(4))
	print("gen tet:",tet.G.generators)
	
	for e in g.G.elements:
		print(p*e)	
