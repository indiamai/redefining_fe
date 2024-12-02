from redefining_fe import *
from test_convert_to_fiat import create_cg1, create_dg1, construct_cg3, construct_rt, construct_nd
import sympy as sp


def test_permute_dg1():
    cell = Point(1, [Point(0), Point(0)], vertex_num=2)

    dg1 = create_dg1(cell)

    for dof in dg1.generate():
        print(dof)

    for g in dg1.cell.group.members():
        print("g", g)
        for dof in dg1.generate():
            print(dof, "->", dof(g))


def test_permute_cg1():
    cell = Point(1, [Point(0), Point(0)], vertex_num=2)

    cg1 = create_cg1(cell)

    for dof in cg1.generate():
        print(dof)

    for g in cg1.cell.group.members():
        print("g", g)
        for dof in cg1.generate():
            print(dof, "->", dof(g))


def test_permute_cg3():
    cell = n_sided_polygon(3)

    cg3 = construct_cg3(cell)

    for dof in cg3.generate():
        print(dof)

    for g in cg3.cell.group.members():
        print(g)
        for dof in cg3.generate():
            print(dof, "->", dof(g))


def test_permute_rt():
    cell = n_sided_polygon(3)

    rt = construct_rt(cell)
    x = sp.Symbol("x")
    y = sp.Symbol("y")
    func = MyTestFunction(sp.Matrix([x, -1/3 + 2*y]), symbols=(x, y))

    for dof in rt.generate():
        print(dof)

    for g in rt.cell.group.members():
        print(g.numeric_rep())
        for dof in rt.generate():
            print(dof, "->", dof(g), "eval, ", dof(g).eval(func))


def test_permute_nd():
    cell = n_sided_polygon(3)

    rt = construct_nd(cell)
    x = sp.Symbol("x")
    y = sp.Symbol("y")
    func = MyTestFunction(sp.Matrix([x, -1/3 + 2*y]), symbols=(x, y))

    for dof in rt.generate():
        print(dof)

    for g in rt.cell.group.members():
        print(g.numeric_rep())
        for dof in rt.generate():
            print(dof, "->", dof(g), "eval, ", dof(g).eval(func))
