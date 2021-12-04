import math
from copy import copy, deepcopy
from print_graphics import print_graphics, plotRadius, plotDistance

import sympy as sp
from sympy import Interval
import numpy as np

def getX(y1, y2, x):
    return sp.solve(sp.Eq(y1, y2), x)[0]

def getY(y1, y2, x):
    ans = sp.solve(sp.Eq(y1, y2), x)
    return y1.subs(x, ans[0])



x = sp.Symbol('x')
y1 = 1.5 - 2 * x
y2 = 2 - 2 * x
y3 = x * 1
y4 = 1.5 * x

print('x1_low', getY(y1, y3, x))
print('x1_high', getY(y2, y4, x))
print('x2_low', getX(y1, y4, x))
print('x2_high', getX(y2, y3, x))

X_old = [[Interval(getX(y1, y4, x), getX(y2, y3, x))], [Interval(getY(y1, y3, x), getY(y2, y4, x))]]
print(X_old)

def sum(i1, i2):
    min = i1.inf + i2.inf
    max = i1.sup + i2.sup
    return Interval(min, max)

def sub(i1, i2):
    min = i1.inf - i2.sup
    max = i1.sup - i2.inf
    return Interval(min, max)

def mul(i1, i2):
    first_data = [i1.inf, i1.sup]
    second_data = [i2.inf, i2.sup]
    muls = []
    for first in first_data:
        for second in second_data:
            muls.append(first * second)
    return Interval(min(muls), max(muls))

def div(i1: Interval, i2: Interval):
    first_data = [i1.inf, i1.sup]
    second_data = [i2.inf, i2.sup]
    divs = []
    for first in first_data:
        for second in second_data:
            divs.append(first / second)
    return Interval(min(divs), max(divs))

def getMid(interval: Interval):
    return (interval.sup + interval.inf) / 2

def getJac(X):
    interval1 = Interval(1, 1)
    interval1rev = Interval(-1, -1)
    interval2 = Interval(2, 2)
    return[[interval1, interval2], [div(interval1, X[1][0]), mul(interval1rev, div(X[0][0], mul(X[1][0], X[1][0])))]]

def getMidMatrix(matrix):
    lam = []
    for i in range(len(matrix)):
        lam.append([])
        for elem in matrix[i]:
            lam[i].append(getMid(elem))
    return np.array(lam, dtype=np.float64)

def getIntervalMatrix(matrix):
    result = []
    for i in range(len(matrix)):
        result.append([])
        for elem in matrix[i]:
            result[i].append(Interval(elem, elem))
    return result

def subMatrix(A, B):
    result = deepcopy(A)
    for i in range(len(result)):
        for j in range(len(result[i])):
            result[i][j] = sub(result[i][j], B[i][j])
    return result

def sumMatrix(A, B):
    result = deepcopy(A)
    for i in range(len(result)):
        for j in range(len(result[i])):
            result[i][j] = sum(result[i][j], B[i][j])
    return result

def mulMatrix(A, B):
    result = []
    for i in range(len(A)):
        result.append([])
        for j in range(len(B[i])):
            result[i].append(Interval(.0, .0))
            for q in range(len(B)):
                result[i][j] = sum(result[i][j], mul(A[i][q], B[q][j]))
    return result

def equal(newX, oldX):
    eps = 10 ** -10
    for i in range(len(newX)):
        if math.fabs(newX[i][0].inf - oldX[i][0].inf) > eps or math.fabs(newX[i][0].sup - oldX[i][0].sup) > eps:
            return False
    return True

def intersection(newX, oldX):
    result = []
    for i in range(len(newX)):
        interval = Interval(max([newX[i][0].inf, oldX[i][0].inf]), min([newX[i][0].sup, oldX[i][0].sup]))
        result.append([interval])
    return result

def Function(X):
    inter1 = sub(sum(X[0][0], mul(Interval(2., 2.), X[1][0])), Interval(1.5, 2.))
    inter2 = sub(div(X[0][0], X[1][0]), Interval(1., 1.5))
    return [[inter1], [inter2]]


b = [[Interval(1.5, 2)], [Interval(1, 1.5)]]
Jac = getJac(X_old)
L = getIntervalMatrix(np.linalg.inv(getMidMatrix(Jac)))
I = [[Interval(1, 1), Interval(0, 0)], [Interval(0, 0), Interval(1, 1)]]
x_mid = getIntervalMatrix(getMidMatrix(X_old))
iterations = []

first_part = subMatrix(x_mid, mulMatrix(L, Function(x_mid))) #x_mid - L * F(x_mid)
second_part = mulMatrix(subMatrix(I, mulMatrix(L, Jac)), subMatrix(X_old, x_mid)) #(I - L * Jac) * (X - x_mid)
K = subMatrix(first_part, second_part) #K(X, x_mid)
X_new = intersection(K, X_old)

while X_new[0][0].is_EmptySet is None and X_new[1][0].is_EmptySet is None and not equal(X_new, X_old):
    Jac = getJac(X_new)
    L = getIntervalMatrix(np.linalg.inv(getMidMatrix(Jac)))
    iterations.append(X_old)
    X_old = deepcopy(X_new)

    x_mid = getIntervalMatrix(getMidMatrix(X_old))
    first_part = subMatrix(x_mid, mulMatrix(L, Function(x_mid)))  # x_mid - L * F(x_mid)
    second_part = mulMatrix(subMatrix(I, mulMatrix(L, Jac)), subMatrix(X_old, x_mid))  # (I - L * Jac) * (X - x_mid)
    K = subMatrix(first_part, second_part)  # K(X, x_mid)
    X_new = intersection(K, X_old)

    print("K: ", K)
    print("X_new: ", X_new)
    print()
print_graphics([])
print_graphics(iterations)
plotRadius(iterations)
plotDistance(iterations)
