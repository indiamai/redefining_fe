

class LinearFunctional():

    def __init__(self, K, A, entity):

        self.K = K
        self.A = A
        self.entity = entity

    def transform(self):
        return LinearFunctional(self.K, self.A.transform(),
                                self.entity.transform())

    def eval(self, v):
        return integrate(self.K * self.A(v), self.entity)


def integrate(exp, entity):
    return None
