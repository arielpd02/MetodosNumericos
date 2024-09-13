
from tools import *



class GaussException(Exception):
    """El sistema no tiene solucion o tiene infinitas soluciones"""


def solve_con_pivot(A: np.array, b: np.array, epsilon: float = 10e-8) -> np.array:
    n, _ = A.shape

    # Agregamos la columna de resultados
    A = np.insert(A, n, b, axis=1)

    # Buscamos el primer pivot
    pivot = 0
    for i in range(n):
        if abs(A[i, 0]) > abs(A[pivot, 0]):
            pivot = i

    # Escalonamos la matriz hacia adelante
    for i in range(n):

        # Pivoteamos si hace falta y seteamos el próximo pivot
        if pivot != i:
            A[[i, pivot]] = A[[pivot, i]]  # O(1)
        pivot = i + 1

        # Si el mayor pivot en valor absoluto es 0, eran todos 0.
        if abs(A[i, i]) == 0: raise GaussException()

        # Si es menor al epsilon de tolerancia, reportamos que hay posible error numérico.
        if abs(A[i, i]) < epsilon: print(f'Posible error numérico:{A[i,i]}s')

        for j in range(i + 1, n):
            # Escalonamos
            A[j, :] = A[j, :] - A[i, :] * (A[j, i] / (A[i, i]))

            # Buscamos el pivot del proximo paso
            if abs(A[j, i + 1]) > abs(A[i + 1, i + 1]):
                pivot = j

    # Despejamos hacia atrás
    for i in range(n - 1, -1, -1):
        for j in range(i + 1, n):
            A[i, n] = A[i, n] - A[i, j] * A[j, n]
        A[i, n] /= A[i, i]

    return A[:, n]
