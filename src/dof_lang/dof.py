
from cell_complex.cells import Point, Edge


class Pairing():

    def __init__(self, entity, space):
        self.entity = entity
        self.space = space


class DeltaPairing(Pairing):

    def __init__(self, entity, space):
        super(DeltaPairing, self).__init__(entity, space)

    def __call__(self, kernel, v):
        assert isinstance(kernel, DeltaKernel)
        print("evaluating v at", kernel.pt)
        print(v)
        return v(*kernel.pt)

    def __repr__(self, kernel):
        return "v(%s)" % str(kernel)
    
    def attach(self, E, attachment):
        return DeltaPairing(E, self.space)


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
        return "integral_{1}({0} * v dx)".format(str(kernel), str(self.entity))
    
    def immerse(self, E):
        raise NotImplementedError("how does the cell transform")


class DeltaKernel():

    def __init__(self, x):
        self.pt = x

    def __repr__(self):
        x = list(map(str, list(self.pt)))
        return ','.join(x)
    
    def immerse(self, attachment, E, V):
        return DeltaKernel(V.pullback(attachment(self.pt)))


class TangentKernel():

    def __init__(self, E):
        self.tangent = E.basis_vectors(return_coords=True)[0]

    def __repr__(self):
        return str(self.tangent)
    
    def trace(self, attachment, E, V):
        return TangentKernel(E)


class DOF():

    def __init__(self, pairing, kernel):
        self.pairing = pairing
        self.kernel = kernel
        self.immersed = False

    def __call__(self, fn):
        if self.immersed:
            return self.pairing(self.kernel,
                                self.pullback(fn.attach(self.attachment)))
        return self.pairing(self.kernel, fn)

    def __repr__(self):
        return self.pairing.__repr__(self.kernel)
    
    def immerse(self, attachment, pullback):
        if not self.immersed:
            self.attachment = attachment
            self.pullback = pullback
        else:
            old_attach = self.attachment
            old_pullback = self.pullback
            self.attachment = lambda x: attachment(old_attach(x))
            self.pullback = lambda v: pullback(old_pullback(v))

        self.immersed = True
    
    # def trace(self, attachment, cell):

    #     return DOF(self.pairing,
    #                self.kernel.trace(attachment, cell))


class MyTestFunction():

    def __init__(self, eq, attach_func = None):
        self.eq = eq
        self.attach_func = None

    def __call__(self, *x):
        print(x)
        print("attach", self.attach_func)
        if self.attach_func:
            print(self.attach_func(tuple(*x)))
            return self.eq(*self.attach_func(tuple(*x)))
        else:
            return self.eq(*x)
    
    def attach(self, attachment):
        print("Attached")
        print(attachment())
        if not self.attach_func:
            new_attach = attachment
        else:
            old_attach = self.attach_func
            new_attach = lambda x: attachment(old_attach(x))
        return MyTestFunction(self.eq, attach_func=new_attach)

    def __repr__(self):
        if self.attach_func:
            return "v(G(x))"
        else:
            return "v(x)"


def construct_point_eval(x, E, V):
    delta_x = DeltaKernel(x)
    return DOF(DeltaPairing(E, V), delta_x)


def construct_tangent_dof(E, V):
    return DOF(L2InnerProd(E, V), TangentKernel(E))
