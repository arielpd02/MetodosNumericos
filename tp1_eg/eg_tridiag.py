import numpy as np
import copy
import matplotlib.pyplot as plt
import sys
from numpy import *
from numpy.linalg import *

#Algortimo de EG para sistemas tridiagonales

def precompute_tridiag(a: np.array, b: np.array, c: np.array) -> np.array:
    n = len(b)
    P = np.zeros((3, n), dtype=float)
    P[1, :] = b
    for i in range(1, n):
        P[1, i] = P[1, i] - (a[i] / P[1, i - 1]) * c[i - 1]
        P[0, i] = a[i] / P[1, i - 1]
    for i in range(n - 2, -1, -1):
        P[2, i] = c[i] / P[1, i + 1]
    return P


def solve_tridiag(a: np.array, b: np.array, c: np.array, d: np.array) -> np.array:
    return solve_tridiag_fast(precompute_tridiag(a, b, c), d)


def solve_tridiag_fast(P: np.array, d: np.array):
    x = np.copy(d)
    for i in range(0, len(x)):
        x[i] = x[i] - P[0, i] * x[i - 1]
    for i in range(len(x) - 2, -1, -1):
        x[i] = x[i] - P[2, i] * x[i + 1]

    return x * (1 / P[1, :])
