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
		for 

if __name__== "__main__":
	g = Group(SymmetricGroup(3))
	print(g.G.is_group)
	print(g.G.conjugacy_classes())
	print(g.G.elements)
	print(g.G.identity)
	g2 = Group(SymmetricGroup(2))
	print(g2.G.is_subgroup(g.G, strict = False))
	p = Permutation([0,2,1])
	print(g.G.contains(p))
	g3 = g.G.subgroup([p])
	print(g.G.coset_factor(p))

	perm = PermutationGroup([p])
	N= Group(perm)
	g.quotient(N)
