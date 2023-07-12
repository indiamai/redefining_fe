from petsc4py import PETSc
import numpy as np

class CellComplex():
	def __init__(self):
		self.group = None
		self.coord_sys = None


class Point(CellComplex):
	
	def __init__(self):
		self.dmplex = PETSc.DMPlex().create()
		self.dmplex.setChart(0,1)
		self.dmplex.setConeSize(0,0)
		super(CellComplex).__init__()


class Interval(CellComplex):
	
	def __init__(self):
		self.dmplex = PETSc.DMPlex().createFromCellList(1, [[0,1]], [[-1.],[1.]])
#		self.dmplex = PETSc.DMPlex().create()
#		self.dmplex.setChart(0, 7)
		#for i in range(2):
		#	self.dmplex.setConeSize(i,0)
		#for i in range(3,6):
	#		self.dmplex.setConeSize(i,2)
	#	self.dmplex.setConeSize(6,3)
	#	self.dmplex.setUp()
	#	self.dmplex.setCone(6, [3,4,5])
	#	self.dmplex.setCone(3, [0,1])
	#	self.dmplex.setCone(4, [0,3])
	#	self.dmplex.setCone(5, [1,2])
		super(CellComplex).__init__()

class Triangle(CellComplex):
	
	def __init__(self):
		self.dmplex = PETSc.DMPlex().createFromCellList(2, [[0,1,2]], [[-1,-np.sqrt(3)/3],[1.,-np.sqrt(3)/3], [0, 2*np.sqrt(3)/3]])
		super(CellComplex).__init__()
stdout = PETSc.Viewer().STDOUT()
p = Point()
i = Interval()
t = Triangle()
print(t.dmplex.getChart())
i.dmplex.topologyView(stdout)

print("Point")
print(p.dmplex.getCone(0))
print("Interval")
for j in range(3):
	print(i.dmplex.getCone(j))
print("Triangle")
for j in range(7):
	print(t.dmplex.getCone(j))
