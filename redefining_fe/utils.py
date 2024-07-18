import numpy as np


def fold_reduce(func_list, *prev):
    """
    Right to left function comprehension

    :param: func_list: list of functions
    :param: prev: starting value(s)
    """
    for func in reversed(func_list):
        prev = func(*prev)
    return prev


def sympy_to_numpy(array, symbols, values):
    """
    Convert a sympy array to a numpy array.

    :param: array: sympy array
    :param: symbols: array of symbols contained in the sympy exprs
    :param: values: array of values to replace the symbols with.

    Due to how sympy handles arrays, we need to squeeze if resulting array
    is greater than 1 dimension to remove extra dimensions
    """
    substituted = array.subs({symbols[i]: values[i] for i in range(len(values))})
    nparray = np.array(substituted).astype(np.float64)

    if len(nparray.shape) > 1:
        return nparray.squeeze()

    if len(nparray.shape) == 0:
        return nparray.item()

    return nparray


def tabulate_sympy(expr, pts):
    res = []
    syms = ["x", "y", "z"]
    for pt in pts:
        if not hasattr(pt, "__iter__"):
            pt = (pt,)
        res.append(expr.subs({syms[i]: pt[i] for i in range(len(pt))}))

    return np.array(res).astype(np.float64)
