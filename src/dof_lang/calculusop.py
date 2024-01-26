class CalculusOp():

    def __init__(self):
        pass

    def transform(self):
        pass

    def evaluate(self, v):
        return "evaluated" + str(v)

    def __call__(self, v):
        return self.evaluate(v)


class HDiv(CalculusOp):
    def __init__(self):
        super(HDiv, self).__init__()

    def evaluate(self, v):
        return "hdiv" + str(v)


if __name__ == "__main__":
    hdiv = HDiv()
    print(hdiv(6))
