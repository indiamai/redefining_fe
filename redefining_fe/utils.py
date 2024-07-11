import numpy as np


# TODO make these the same function
def fold_reduce_group(func_list, x):
    """ nested function composition helper function, right to left """
    prev = x
    for func in reversed(func_list):
        prev = func(prev)
    return prev


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
    substituted = array.subs({symbols[i]: values[i] for i in range(len(values))})
    nparray = np.array(substituted).astype(np.float64)

    if len(nparray.shape) > 1:
        return nparray.squeeze()
    return nparray
