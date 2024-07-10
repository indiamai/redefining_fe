from redefining_fe import *
# .. plot:: ../../test/plots.py


def construct_dg1():
    # [test_dg1_int 0]
    edge = Point(1, [Point(0), Point(0)], vertex_num=2)
    xs = [DOF(DeltaPairing(), PointKernel((-1,)))]
    dg1 = ElementTriple(edge, (P1, CellL2, C0), DOFGenerator(xs, S2, S1))
    # [test_dg1_int 1]
    return dg1


def plot_dg1():
    dg1 = construct_dg1()
    dg1.plot()


plot_dg1()
