import math
from copy import deepcopy

A0 = [[0.20, 1.60, -0.11],
      [0.20, -0.10, 0.90],
      [-0.50, -0.20, -0.31]]
b0 = [1.98, 2.30, -2.32]


def Gauss_method(A, b):
    for i in range(len(A)):
        a = A[i][i]
        for k in range(i, len(A)):
            A[i][k] /= a
        b[i] /= a
        for k in range(i + 1, len(A)):
            c = -A[k][i]
            for j in range(i, len(A)):
                A[k][j] += c * A[i][j]
            b[k] += c * b[i]

    for i in range(len(A) - 1, -1, -1):
        for k in range(i - 1, -1, -1):
            c = -A[k][i]
            A[k][i] += c * A[i][i]
            b[k] += c * b[i]

    return b.copy()


def Gauss_main_value(A, b):
    n = len(A)
    visited = [False] * n
    order = []
    for _ in range(n):
        max_items = map(lambda x: (x[0], *max(enumerate(x[1]), key=lambda y: abs(y[1]))), enumerate(A))
        i, j = max(filter(lambda y: not visited[y[0]], max_items), key=lambda x: abs(x[2]))[:2]
        a = A[i][j]
        visited[i] = True
        order.append((i, j))
        for k in range(n):
            A[i][k] /= a
        b[i] /= a
        for k in range(n):
            if visited[k]:
                continue
            c = -A[k][j]
            for l in range(n):
                A[k][l] += c * A[i][l]
            b[k] += c * b[i]

    answer = [None] * n
    for i in range(n - 1, -1, -1):
        for k in range(i - 1, -1, -1):
            c = -A[order[k][0]][order[i][1]]
            A[order[k][0]][order[i][1]] += c * A[order[i][0]][order[i][1]]
            b[order[k][0]] += c * b[order[i][0]]
        answer[order[i][1]] = b[order[i][0]]

    return answer


def get_distance(x, y):
    return math.sqrt(sum([(x[i] - y[i]) ** 2 for i in range(len(x))]))


def main():
    A = deepcopy(A0)
    b = deepcopy(b0)
    res = [3.0, 1.0, 2.0]

    res1 = Gauss_method(A, b)
    A = deepcopy(A0)
    b = deepcopy(b0)
    res2 = Gauss_main_value(A, b)
    print("Точное решение ", res)
    print("Решение методом Гаусса:")
    print(res1)
    print("Погрешность: ", get_distance(res, res1))
    print("-"*60)
    print("Решение методом Гаусса c выбором главного элемента:")
    print(res2)
    print("Погрешность: ", get_distance(res, res2))


if __name__ == "__main__":
    main()
