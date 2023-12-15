from petsc4py import PETSc
import numpy as np
from enum import Enum
from sympy.combinatorics import Permutation
from groups.group import S1,S2,S3,S4 
import matplotlib.pyplot as plt

class CellComplex():
	"""
	Abstract Class for cell complexes

	Defines operational structure of cell complexes

	* TODO: make this actually an abstract class?
	"""
	def __init__(self):
		self.points = {}
		self.edges ={} 
		self.dmplex = None

	def add_point(self, point):
		assert isinstance(point, Point)
		num = len(self.points) 
		self.points[num] = point

	
	def connect(self, lower, upper, immersion, orientation=None):
		assert self.points[lower].dim < self.points[upper].dim
		if not orientation:
			orientation = self.points[lower].group.G.identity
		assert self.points[lower].group.is_element(orientation)
		self.edges[upper, lower] = (immersion, orientation)	
	
	def construct_dmplex(self):
		self.dmplex = PETSc.DMPlex().create()
		self.dmplex.setChart(0, len(self.points))		
		cones = {}
		for u,l in self.edges:
			if u in cones.keys():
				cones[u].append(l)
			else:
				cones[u] = [l]	
		
		for i in cones.keys():
			self.dmplex.setConeSize(i, len(cones[i]))	
		self.dmplex.setUp()
		for i in cones.keys():
		        self.dmplex.setCone(i, cones[i])
			

	def get_immersion_operator(self,initial_level, output_level, entity_association):
		# note entity association is currently given by the numbering in the dmplex but may need to look at this
		# function could have an added option to check all possible paths? that might be too complex. need to assert all paths should be the same.
		# this should also handle orientations? 
		immersion_ops = []
		current_entity = entity_association
		for i in range(output_level - initial_level):
			ops = [(upper, self.edges[upper, current_entity][0]) for upper in range(len(self.points)) if (upper, current_entity) in self.edges]
			current_entity = ops[0][0]
			immersion_ops.append(ops[0][1])

		#hacky way to nest the functions
		def immersion_op(x):
			res = x
			for i in range(len(immersion_ops), 0, -1):
				res = immersion_ops[i-1](res)
			return res
		return immersion_op

	def plot(self):
		counter = {}
		coords = {}
		for p in self.points:
			point = self.points[p]
			if point.dim not in counter.keys():
				counter[point.dim] = 0
			coords[p] = (counter[point.dim], point.dim)
			plt.plot(counter[point.dim], point.dim,  "bo", label=p)
			plt.annotate(p, coords[p], fontsize=20)
			counter[point.dim] += 1	
		
		for (u,l) in self.edges:
			(x1,y1) = coords[u]
			(x2,y2) = coords[l]
			plt.plot([x1,x2], [y1,y2], color="black")
			plt.annotate(self.points[l].group.perm_name(self.edges[u,l][1]), ((x1+x2)/2, (y1+y2)/2))
			
		plt.axis('off')		
		plt.show()

class Point():
	"""
	Class to represent Point types
	
	"""
	dim = 0
	def __init__(self, group):
		self.group = group		

class Vertex(Point):
	"""
	Dimension: 0
	Coordinate System: Empty (no points)
	"""

	def __init__(self):
		self.dim = 0
		self.coord = ()
		super(Vertex, self).__init__(S1())


class Interval(Point):
	"""
	Class to represent interval points
	Dimension 1
	Coordinate System: (x)
	"""

	def __init__(self):
		self.dim = 1	
		super(Interval, self).__init__(S2())

	def default_cell_complex(self):
		cell = CellComplex()
		cell.add_point(Interval())
		cell.add_point(Vertex())
		cell.add_point(Vertex())
		i_0 = lambda v: v(-1)
		i_1 = lambda v: v(1)
		cell.connect(2, 0, i_0)
		cell.connect(1, 0, i_1)
		return cell
		
	

class Triangle(Point):
	""" 
	Class to represent triangle points 
	Dimension 2
	Coordinate System: (x,y)
	"""
	def __init__(self):
		self.dim = 2
		super(Triangle, self).__init__(S3())

	def default_cell_complex(self):
		cell = CellComplex()
		cell.add_point(Triangle())
		cell.add_point(Vertex())
		cell.add_point(Vertex())
		cell.add_point(Vertex())
		cell.add_point(Interval())
		cell.add_point(Interval())
		cell.add_point(Interval())
		r = Permutation([1,0])
		i_v0 = lambda v: v(-1)
		i_v1 = lambda v: v(1)
		i_e0 = lambda v: (lambda x: v(x, -np.sqrt(3)/3)) 
		i_e1 = lambda v: (lambda a: v(-1+(a + 1)/2, -np.sqrt(3)/3 + ((a+1)/2)*np.sqrt(3)))
		i_e2 = lambda v: (lambda a: v(1-(a+1)/2, -np.sqrt(3)/3 + ((a+1)/2)*np.sqrt(3)))
		cell.connect(4,0, i_e0)
		cell.connect(5,0, i_e1)
		cell.connect(6,0, i_e2)
		cell.connect(1,4, i_v0)
		cell.connect(2,4, i_v1)
		cell.connect(2,5, i_v0)
		cell.connect(3,5, i_v1)
		cell.connect(1,6, i_v0)
		cell.connect(3,6, i_v1)
		return cell

	def immersion_e0(self, v):
		# v is a function with one argument
		# output is a function of two arguments (???)
                pass

class Tetrahedron(Point):
	""" 
	Class to represent triangle points 
	Dimension 3
	Coordinate System: (x,y,z)
	TODO: Check symmetry group of tets
	"""
	def __init__(self):
		self.dim = 3
		super(Tetrahedron, self).__init__(S4())



if __name__== "__main__":
	m = Interval()
	e = m.group.G.identity	
	a = Vertex()
	cell = m.default_cell_complex() 
	cell.construct_dmplex()
#	cell.plot()
	imm = cell.get_immersion_operator(0,1,1)
	print("entity 1")
	print(imm(lambda x: x**2 + 1))
	imm = cell.get_immersion_operator(0,1,2)
	print("entity 2")
	print("output", imm(lambda x: x**2 + 1))
	cell2 =Triangle().default_cell_complex() 
#	cell2.plot()
	cell2.construct_dmplex()
	i1 = cell2.get_immersion_operator(0,1,2)
	i2 = cell2.get_immersion_operator(0,2,3)
	i4 = cell2.get_immersion_operator(1,2,4)
	i5 = cell2.get_immersion_operator(1,2,5)
	i6 = cell2.get_immersion_operator(1,2,6)
	print(i1(lambda x: x))
	print(i2(lambda x,y:(x,y)))
	space = np.linspace(-1,1,20)
	edge0 = ([i4(lambda x,y: (x,y))(s) for s in space])
	unzipped0 = list(zip(*edge0))
	plt.plot(unzipped0[0], unzipped0[:][1])
	edge1 = [i6(lambda x,y: (x,y))(s) for s in space]
	unzipped1 = list(zip(*edge1))
	plt.plot(unzipped1[0], unzipped1[1])
	edge2 = [i5(lambda x,y: (x,y))(s) for s in space]
	unzipped2 = list(zip(*edge2))
	plt.plot(unzipped2[0], unzipped2[1])
	plt.show()
	cell3 = CellComplex()
	cell3.add_point(Tetrahedron())
	cell3.add_point(Triangle())
	cell3.add_point(Triangle())
	cell3.add_point(Triangle())
	print(cell3.points)
	rot = Permutation([2,0,1])
	r = Permutation([0,2,1])
	cell3.connect(1,0, "i3", rot)
	cell3.connect(2,0, "i3", r)
	cell3.connect(3,0, "i3")
#	cell3.get_immersion_operator(0,2,2)
#	cell3.plot()
