
class ElementSpaceTuple():

    def __init__(self, v, sobolev, interp_domain):
        assert isinstance(v, ContinuousFunctionSpace)
        assert isinstance(sobolev, ContinuousFunctionSpace)
        assert isinstance(interp_domain, ContinuousFunctionSpace)

        self.v = v
        self.w = sobolev
        self.w_I = interp_domain


class ContinuousFunctionSpace():

    def __init__(self, name, domain=None):
        self.name = name

    def trace(self, v):
        raise NotImplementedError("Trace not implemented for", str(type(self)))
    

class H1(ContinuousFunctionSpace):

    def __init__(self):
        super(H1, self).__init__("H1")

    def pullback(self, v):
        # temporarily everything is reference spacE?
        return v


class L2(ContinuousFunctionSpace):

    def __init__(self):
        super(L2, self).__init__("L2")

    def pullback(self, v):
        # temporarily everything is reference spacE?
        return v

if __name__ == "__main__":
    h1 = ContinuousFunctionSpace("H1")
    h1.trace()
