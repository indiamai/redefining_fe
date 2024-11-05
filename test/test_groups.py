from pytest import *
from redefining_fe import *


def test_numerical_orientations():
    vert = Point(0)
    print(vert.group.compute_num_reps())
    edge = Point(1, [Point(0), Point(0)], vertex_num=2)
    # group = S2.add_cell(edge)
    print(edge.group.compute_num_reps())

    cell = n_sided_polygon(3)
    # group = S3.add_cell(cell)
    print(cell.group.compute_num_reps())
    mems = cell.group.members()
    for m in mems:
        print(m.compute_num_rep())


# def test_group_equality():
#     cell = n_sided_polygon(3)

#     s1 = S1.add_cell(cell)
#     s1_new = S1.add_cell(cell)

#     assert s1 == s1_new

    # cell2 = Point(1, [Point(0), Point(0)], vertex_num=2)

    # s1 = S2.add_cell(cell)
    # s1_new = S2.add_cell(cell2)

    # assert not s1 == s1_new
