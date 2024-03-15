
from cell_complex.cells import Point, Edge
from FIAT.quadrature import GaussLegendreQuadratureLineRule

from FIAT.reference_element import DefaultLine

import numpy as np

class Pairing():

    def __init__(self, entity, space):
        self.entity = entity
        self.space = space


class DeltaPairing(Pairing):

    def __init__(self, entity, space):
        super(DeltaPairing, self).__init__(entity, space)

    def __call__(self, kernel, v):
        assert isinstance(kernel, PointKernel)
        return v(*kernel.pt)

    def __repr__(self, kernel, fn):
        return "{0}({1})".format(fn, str(kernel))


class L2InnerProd(Pairing):
    """ need to think about the abstraction level here - 
    are we wanting to define them as quadrature now? or defer this?
    """
    def __init__(self, entity, space):
        super(L2InnerProd, self).__init__(entity, space)

    def __call__(self, kernel, v):
        # evaluates integral
        print("evaluating", kernel, v)
        assert self.entity.dim() == 1
        quadrature = GaussLegendreQuadratureLineRule(DefaultLine(), 5)
        # print(kernel())
        # print(np.dot(kernel(), v(-0.9061)))
        return quadrature.integrate(lambda *x: np.dot(kernel(), v(*x)), unpack=True)

    def __repr__(self, kernel, fn):
        return "integral_{1}({0} * {2}) dx)".format(str(kernel),
                                                    str(self.entity),
                                                    fn)


class PointKernel():

    def __init__(self, x):
        self.pt = x

    def __repr__(self):
        x = list(map(str, list(self.pt)))
        return ','.join(x)
    
    def __call__(self):
        return self.pt


class GeneralKernel():

    def __init__(self, fn):
        self.fn = fn

    def __repr__(self):
        return str(self.fn)

    def __call__(self):
        return self.fn


class DOF():

    def __init__(self, pairing, kernel):
        self.pairing = pairing
        self.kernel = kernel
        self.immersed = False

    def __call__(self, fn):
        if self.immersed:
            attached_fn = fn.attach(self.attachment)
            return self.pairing(self.kernel,
                                self.target_space.pullback(attached_fn, self.trace_entity))
        return self.pairing(self.kernel, fn)

    def __repr__(self):
        if self.immersed:
            fn = "tr_{1}_{0}(v)".format(str(self.trace_entity),
                                        str(self.target_space))
        else:
            fn = "v"
        return self.pairing.__repr__(self.kernel, fn)

    def immerse(self, entity, attachment, target_space):
        if not self.immersed:
            self.trace_entity = entity
            self.attachment = attachment
            self.target_space = target_space
        else:
            raise RuntimeError("Error immersing twice")
            # old_attach = self.attachment
            # old_pullback = self.pullback
            # self.trace_entity = entity
            # self.attachment = lambda x: attachment(old_attach(x))
            # self.pullback = lambda v: pullback(old_pullback(v))

        self.immersed = True


class MyTestFunction():

    def __init__(self, eq, attach_func=None):
        self.eq = eq
        self.attach_func = attach_func

    def __call__(self, *x):
        if self.attach_func:
            return self.eq(*self.attach_func(*x))
        else:
            return self.eq(*x)

    def attach(self, attachment):
        if not self.attach_func:
            return MyTestFunction(self.eq, attach_func=attachment)
        else:
            old_attach = self.attach_func
            return MyTestFunction(self.eq,
                                  attach_func=lambda *x: attachment(old_attach(*x)))

    def __repr__(self):
        if self.attach_func:
            return "v(G(x))"
        else:
            return "v(x)"
