

class Generator():

    def __init__(self, functional):

        self.K = K
        self.A = A
        self.entity = entity
    
    def __call__(self, g):
        

    def transform(self):
        return LinearFunctional(self.K, self.A.transform(),
                                self.entity.transform())

    def eval(self, v):
        return integrate(self.K * self.A(v), self.entity)


def integrate(exp, entity):
    return None
