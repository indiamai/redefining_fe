from petsc4py import PETSc
import numpy as np
from enum import Enum

class Group(Enum):
	"""
	Enum to temporarily represent the symmetry groups
	"""
	S1 = 1
	S2 = 2
	C3 = 3
	S3 = 4
	D4 = 5
	S4 = 6


class CellComplex():
	"""
	Abstract* Class for cell complexes

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
	
	def connect(self, upper, lower, immersion, orientation):
		assert self.points[upper].dim < self.points[lower].dim
		self.edges[lower, upper] = (immersion, orientation)	
	
	def construct_dmplex(self):
		self.dmplex = PETSc.DMPlex().create()
		self.dmplex.setChart(0, len(self.points))		
		cones = {}
		for u,l in self.edges:
			if u in cones.keys():
				cones[u].append(l)
			else:
				cones[u] = [l]	
		
		for cone in cones:
			print(cone)		
		
		
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
		super(Vertex, self).__init__(Group.S1)

class Interval(Point):
	"""
	Class to represent interval points
	Dimension 1
	Coordinate System: (x)
	"""

	def __init__(self):
		self.dim = 1
		super(Interval, self).__init__(Group.S2)
	

class Triangle(CellComplex):
	""" 
	Class to represent triangle points 
	Dimension 2
	Coordinate System: (x,y)
	"""
	def __init__(self):
		super(Triangle, self).__init__(Group.S3)




if __name__== "__main__":
	a = Vertex()
	print(a.group)
	cell = CellComplex()
	cell.add_point(a)
	cell.add_point(Vertex())
	cell.add_point(Interval())
	cell.connect(0, 2, "immersion", "orientation")
	cell.connect(1, 2, "immersion", "orientation")
	print(cell.edges)
	cell.construct_dmplex()
