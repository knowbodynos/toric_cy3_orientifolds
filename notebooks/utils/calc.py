import numpy as np
from itertools import combinations


def get_srideal(triang):
    triang = np.array(triang)
    n = triang.shape[1]
    k = triang.max() + 1
    srideal = []
    for i in range(n):
        for c in combinations(range(k), i + 1):
            issubset = False
            for sc in srideal:
                if set(sc).issubset(set(c)):
                    issubset = True
                    break
            if issubset:
                continue
            issubset = False
            for t in triang:
                if set(c).issubset(set(t)):
                    issubset = True
                    break
            if not issubset:
                srideal.append(list(c))
    return srideal