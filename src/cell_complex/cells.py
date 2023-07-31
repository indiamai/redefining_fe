from petsc4py import PETSc
import numpy as np
from enum import Enum
from sympy.combinatorics import Permutation
from src.groups.group import S1,S2,S3,S4 
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
			print(self.edges[u,l])
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
		cell.connect(2, 0, "immersion")
		cell.connect(1, 0, "immersion")
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
		cell.connect(4,0, "i2")
		cell.connect(5,0, "i2",r)
		cell.connect(6,0, "i2",r)
		cell.connect(1,4, "i1")
		cell.connect(2,4, "i1")
		cell.connect(2,5, "i1")
		cell.connect(3,5, "i1")
		cell.connect(1,6, "i1")
		cell.connect(3,6, "i1")
		return cell

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
	cell.plot()
	for i in range(3):
		print(i, " ",cell.dmplex.getCone(i))
	cell2 =Triangle().default_cell_complex() 
	cell2.plot()
	cell2.construct_dmplex()
	for i in range(7):
		print(i, " ", cell2.dmplex.getCone(i))

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
	cell3.plot()
