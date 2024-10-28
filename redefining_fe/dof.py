from FIAT.quadrature_schemes import create_quadrature
from FIAT.functional import PointEvaluation
from redefining_fe.cells import CellComplexToFiat
import numpy as np


class Pairing():
    """
    Akin to an inner product, the pairing combines a kernel and an input function
    """

    def __init__(self):
        self.entity = None

    def add_entity(self, entity):
        self.entity = entity

    def _to_dict(self):
        o_dict = {"entity": self.entity}
        return o_dict


class DeltaPairing(Pairing):
    """
    The delta pairing allows the evaluation at a single points

    Calling method:
    :param: kernel: Normally a PointKernel
    """

    def __init__(self):
        super(DeltaPairing, self).__init__()

    def __call__(self, kernel, v):
        assert isinstance(kernel, PointKernel)
        return v(*kernel.pt)

    def convert_to_fiat(self, ref_el, pt):
        # assert isinstance(kernel, PointKernel)
        return PointEvaluation(ref_el, pt)

    def __repr__(self):
        return "{fn}({kernel})"

    def dict_id(self):
        return "Delta"

    def _from_dict(obj_dict):
        new_obj = DeltaPairing()
        new_obj.add_entity(obj_dict["entity"])
        return new_obj


class L2InnerProd(Pairing):
    """ need to think about the abstraction level here -
    are we wanting to define them as quadrature now? or defer this?
    """
    def __init__(self):
        super(L2InnerProd, self).__init__()

    def __call__(self, kernel, v):
        # print("evaluating", kernel, v, "on", self.entity)
        quadrature = create_quadrature(CellComplexToFiat(self.entity), 5)

        def kernel_dot(x):
            return np.dot(kernel(*x), v(*x))
        return quadrature.integrate(kernel_dot)

    def convert_to_fiat(self, ref_el, kernel):
        raise NotImplementedError("L2 functionals not yet fiat convertible")

    def __repr__(self):
        return "integral_{}({{kernel}} * {{fn}}) dx)".format(str(self.entity))

    def dict_id(self):
        return "L2Inner"

    def _from_dict(obj_dict):
        new_obj = L2InnerProd()
        new_obj.add_entity(obj_dict["entity"])
        return new_obj


class PointKernel():

    def __init__(self, x):
        if not isinstance(x, tuple):
            x = (x,)
        self.pt = x

    def __repr__(self):
        x = list(map(str, list(self.pt)))
        return ','.join(x)

    def permute(self, g):
        return PointKernel(g(self.pt))

    def __call__(self, *args):
        return self.pt

    def _to_dict(self):
        o_dict = {"pt": self.pt}
        return o_dict

    def dict_id(self):
        return "PointKernel"

    def _from_dict(obj_dict):
        return PointKernel(tuple(obj_dict["pt"]))


class PolynomialKernel():

    def __init__(self, fn):
        self.fn = fn

    def __repr__(self):
        return str(self.fn)

    def permute(self, g):
        def permuting(*x):
            return self.fn(*g(x))
        return PolynomialKernel(permuting)

    def __call__(self, *args):
        res = self.fn(*args)
        return res

    def _to_dict(self):
        o_dict = {"fn": self.fn}
        return o_dict

    def dict_id(self):
        return "PolynomialKernel"

    def _from_dict(obj_dict):
        return PolynomialKernel(obj_dict["fn"])


class DOF():

    def __init__(self, pairing, kernel, entity=None, attachment=None, target_space=None, g=None, immersed=False, generation=dict({}), sub_id=None):
        self.pairing = pairing
        self.kernel = kernel
        self.immersed = immersed
        self.trace_entity = entity
        self.attachment = attachment
        self.target_space = target_space
        self.g = g
        self.id = None
        self.sub_id = sub_id
        self.generation = generation
        if entity is not None:
            self.pairing.add_entity(entity)

    def __call__(self, g):
        new_generation = self.generation.copy()
        return DOF(self.pairing, self.kernel.permute(g), self.trace_entity, self.attachment, self.target_space, g, self.immersed, new_generation, self.sub_id)

    def eval(self, fn, pullback=True):
        return self.pairing(self.kernel, fn)

    def add_context(self, cell, space, g, overall_id=None, generator_id=None):
        # We only want to store the first instance of each
        if self.trace_entity is None:
            self.trace_entity = cell
            self.generation[self.trace_entity.dim()] = g
            self.pairing.add_entity(cell)
        if self.target_space is None:
            self.target_space = space
        if self.id is None and overall_id is not None:
            self.id = overall_id
        if self.sub_id is None and generator_id is not None:
            self.sub_id = generator_id

    def convert_to_fiat(self, ref_el):
        if isinstance(self.kernel, PointKernel):
            pt = self.eval(MyTestFunction(lambda *x: x))
            return self.pairing.convert_to_fiat(ref_el, pt)
        raise NotImplementedError("Fiat conversion only implemented for Point eval")

    def __repr__(self, fn="v"):
        return str(self.pairing).format(fn=fn, kernel=self.kernel)

    def immerse(self, entity, attachment, target_space, g, triple):
        new_generation = self.generation.copy()
        new_generation[target_space.domain.dim()] = g
        return ImmersedDOF(self.pairing, self.kernel, entity, attachment, target_space, g, triple, new_generation, self.sub_id)

    def _to_dict(self):
        """ almost certainly needs more things"""
        o_dict = {"pairing": self.pairing, "kernel": self.kernel}
        return o_dict

    def dict_id(self):
        return "DOF"

    def _from_dict(obj_dict):
        return DOF(obj_dict["pairing"], obj_dict["kernel"])


class ImmersedDOF(DOF):
    def __init__(self, pairing, kernel, entity=None, attachment=None, target_space=None, g=None, triple=None, generation={}, sub_id=None):
        self.immersed = True
        self.triple = triple
        super(ImmersedDOF, self).__init__(pairing, kernel, entity=entity, attachment=attachment, target_space=target_space, g=g, immersed=True, generation=generation, sub_id=sub_id)

    def eval(self, fn, pullback=True):
        attached_fn = fn.attach(self.attachment)

        if not pullback:
            return self.pairing(self.kernel, attached_fn)

        return self.pairing(self.kernel,
                            self.target_space(attached_fn, self.trace_entity, self.g))

    def __call__(self, g):
        return ImmersedDOF(self.pairing, self.kernel.permute(g), self.trace_entity,
                           self.attachment, self.target_space, g, self.immersed, self.sub_id)

    def __repr__(self):
        fn = "tr_{1}_{0}(v)".format(str(self.trace_entity), str(self.target_space))
        return super(ImmersedDOF, self).__repr__(fn)

    def immerse(self, entity, attachment, trace, g):
        raise RuntimeError("Error: Immersing twice not supported")


class MyTestFunction():

    def __init__(self, eq, attach_func=None, symbols=None):
        self.eq = eq
        self.attach_func = attach_func
        self.symbols = symbols

    def __call__(self, *x, sym=False):
        if self.symbols:
            if self.attach_func and not sym:
                res = self.eq.subs({symb: val for (symb, val) in zip(self.symbols, self.attach_func(*x))})
            else:
                res = self.eq.subs({symb: val for (symb, val) in zip(self.symbols, x)})
            if res.free_symbols == set():
                array = np.array(res).astype(np.float64)
                return array
            else:
                return res
        if self.attach_func and not sym:
            return self.eq(*self.attach_func(*x))
        else:
            # TODO remove this as will already be symbolic
            return self.eq(*x)

    def attach(self, attachment):
        if not self.attach_func:
            return MyTestFunction(self.eq, attach_func=attachment, symbols=self.symbols)
        else:
            old_attach = self.attach_func
            if self.symbols:
                return MyTestFunction(self.eq,
                                      attach_func=attachment(old_attach(*self.symbols)),
                                      symbols=self.symbols)
            else:
                return MyTestFunction(self.eq,
                                      attach_func=lambda *x: attachment(old_attach(*x)))

    def __repr__(self):
        if self.attach_func:
            return "v(G(x))"
        else:
            return "v(x)"

    def _to_dict(self):
        return {"eq": self.eq}

    def dict_id(self):
        return "Function"
