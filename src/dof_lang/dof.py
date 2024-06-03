
from cell_complex.cells import Point, Edge
from FIAT.quadrature import GaussLegendreQuadratureLineRule
from FIAT.reference_element import DefaultLine

import numpy as np

class Pairing():

    def __init__(self):
        self.entity = None
        self.space = None

    def add_entity(self, entity):
        self.entity = entity


class DeltaPairing(Pairing):

    def __init__(self):
        super(DeltaPairing, self).__init__()

    def __call__(self, kernel, v):
        assert isinstance(kernel, PointKernel)
        return v(*kernel.pt)

    def __repr__(self):
        return "{fn}({kernel})"


class L2InnerProd(Pairing):
    """ need to think about the abstraction level here - 
    are we wanting to define them as quadrature now? or defer this?
    """
    def __init__(self):
        super(L2InnerProd, self).__init__()

    def __call__(self, kernel, v):
        # evaluates integral on generic edge only
        # print("evaluating", kernel, v, "on", self.entity)
        if self.entity.dim() == 1:
            quadrature = GaussLegendreQuadratureLineRule(DefaultLine(), 5)
            return quadrature.integrate(lambda *x: np.dot(kernel(), v(*x)), unpack=True)
        elif self.entity.dim() == 2:
            # ("evaluating at origin")
            return np.dot(kernel(), v(0, 0))
        else:
            raise NotImplementedError("L2 Inner product evaluation not implemented for this dimension")

    def __repr__(self):
        return "integral_{}({{kernel}} * {{fn}}) dx)".format(str(self.entity))


class PointKernel():

    def __init__(self, x):
        self.pt = x

    def __repr__(self):
        x = list(map(str, list(self.pt)))
        return ','.join(x)
    
    def permute(self, g):
        return PointKernel(g(self.pt))

    def __call__(self):
        return self.pt


class GeneralKernel():

    def __init__(self, fn):
        self.fn = fn

    def __repr__(self):
        return str(self.fn)

    def permute(self, g):
        return self

    def __call__(self):
        return self.fn


class DOF():

    def __init__(self, pairing, kernel, immersed=False,
                 entity=None, attachment=None, target_space=None, g = None):
        self.pairing = pairing
        self.kernel = kernel
        self.immersed = immersed
        self.trace_entity = entity
        self.attachment = attachment
        self.target_space = target_space
        self.g = g
        if entity is not None:
            self.pairing.add_entity(entity)

    def __call__(self, g):
        return DOF(self.pairing, self.kernel.permute(g), self.immersed,
                   self.trace_entity, self.attachment, self.target_space, g)

    def eval(self, fn, pullback=True):
        if self.immersed:
            attached_fn = fn.attach(self.attachment)

            if not pullback:
                return self.pairing(self.kernel, attached_fn)
            
            return self.pairing(self.kernel,
                                self.target_space.pullback(attached_fn,
                                                           self.trace_entity,
                                                           self.g))
        return self.pairing(self.kernel, fn)
    
    # def fn_eval(self, fn, pullback=True):
    #     if self.immersed:
    #         attached_fn = fn.attach(self.attachment)
    #         if pullback:
    #             return self.target_space.pullback(attached_fn, self.trace_entity,
    #                                                            self.g)
    #         else:
    #             return self.pairing(self.kernel, attached_fn)
    #     return fn

    def add_context(self, cell, space):
        # We only want to store the first instance of each
        if self.trace_entity is None:
            self.trace_entity = cell
            self.pairing.add_entity(cell)
        if self.target_space is None:
            self.target_space = space

    def __repr__(self):
        if self.immersed:
            fn = "tr_{1}_{0}(v)".format(str(self.trace_entity),
                                        str(self.target_space))
        else:
            fn = "v"
        return str(self.pairing).format(fn=fn, kernel=self.kernel)

    def immerse(self, entity, attachment, target_space, g):
        if not self.immersed:
            return DOF(self.pairing, self.kernel,
                       True, entity, attachment, target_space, g)
        else:
            raise RuntimeError("Error: Immersing twice not supported")


class MyTestFunction():

    def __init__(self, eq, attach_func=None):
        self.eq = eq
        self.attach_func = attach_func

    def __call__(self, *x, sym=False):
        if self.attach_func and not sym:
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
