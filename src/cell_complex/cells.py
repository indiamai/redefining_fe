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

	def add_point(self, point, num):
		assert isinstance(point, Point)
		if num not in self.points.keys(): 
			self.points[num] = point
		else:
			print("Point already exists")
	
	def connect(self, upper, lower, immersion, orientation):
		assert self.points[upper].dim < self.points[lower].dim
		self.edges[lower, upper] = (immersion, orientation)	
		
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
	Dimension 1
	Coordinate System: [-1,1]
	"""

	def __init__(self):
		self.dim = 1
		super(Interval, self).__init__(Group.S2)

class Triangle(CellComplex):
	""" Class to represent triangle cells

	Dimension 2

	Hasse Diagram:

		1     2     3
                | / \   /\  | 
		4     6     5 
		 \    |     /
		      0
	"""
	def __init__(self):
		self.dmplex = PETSc.DMPlex().createFromCellList(2, [[0,1,2]], [[-1,-np.sqrt(3)/3],[1.,-np.sqrt(3)/3], [0, 2*np.sqrt(3)/3]])
		super(CellComplex).__init__()




if __name__== "__main__":
	a = Vertex()
	print(a.group)
	cell = CellComplex()
	cell.add_point(a, 1)
	cell.add_point(Vertex(), 2)
	cell.add_point(Interval(), 0)
	cell.connect(1, 0, "immersion", "orientation")
	print(cell.edges)
