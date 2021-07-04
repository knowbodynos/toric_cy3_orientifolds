import numpy as np
from sympy.combinatorics import Permutation


def mat(mat_):
    return np.array(eval(mat_.replace('{', '[').replace('}', ']')))

def invol(invol_, ndivs):
    pairs = set(frozenset(int(y.lstrip('D')) - 1 for y in x.split('->')) for x in invol_.strip('{}').split(','))
    perm = Permutation()
    for x in pairs:
        perm *= Permutation(*x)
    for x in set(range(ndivs)) - set.union(pairs):
        perm *= Permutation(x)
    return perm
