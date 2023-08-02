



class LinearFunctional():

	def __init__(self, K, A, entity):
		
		self.K = K
		self.A = A
		self.entity = entity
	
	def transform(self):
		return LinearFunctional(K, self.A.transform(), entity.transform())
	
	def eval(self, v):
		return integrate(K*A(v), entity)


def integrate(exp, entity):
	return None
