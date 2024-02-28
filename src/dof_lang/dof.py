

from typing import Any


class Pairing():

    def __init__(self, entity, space):
        self.entity = entity
        self.space = space


class DeltaPairing(Pairing):

    def __init__(self, entity, space):
        super(DeltaPairing, self).__init__(entity, space)

    def __call__(self, kernel, v):
        assert isinstance(kernel, DeltaKernel)
        return v(kernel.pt)

    def __repr__(self, kernel):
        return "v(%s)" % str(kernel)


class L2InnerProd(Pairing):
    """ need to think about the abstraction level here - 
    are we wanting to define them as quadrature now? or defer this?
    """
    def __init__(self, entity, space):
        super(L2InnerProd, self).__init__(entity, space)

    def __call__(self, x, v):
        # evaluates integral
        pass

    def __repr__(self, kernel):
        return "integral({0} * v dx)".format(str(kernel))


class DeltaKernel():

    def __init__(self, x):
        self.pt = x

    def __repr__(self):
        x = list(map(str, list(self.pt)))
        return ','.join(x)


class DOF():

    def __init__(self, pairing, kernel):
        self.pairing = pairing
        self.kernel = kernel

    def __call__(self, fn):
        return self.pairing(self.kernel, fn)
    
    def __repr__(self):
        return self.pairing.__repr__(self.kernel)


def construct_point_eval(x, E, V):
    delta_x = DeltaKernel(x)
    return DOF(DeltaPairing(E, V[0]), delta_x)

if __name__ == "__main__":
    delta_0 = DeltaKernel(0)
    print(delta_0)

    eval_0_dof = DOF(DeltaPairing("edge", "h1"), delta_0)
    print(eval_0_dof)

    int_e_dof = DOF(L2InnerProd("edge", "h2"), "x")
    print(int_e_dof)