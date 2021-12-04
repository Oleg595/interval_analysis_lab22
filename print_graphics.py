import math

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as pth

def print_graphics(iterations):
    fig, ax = plt.subplots()
    x = np.linspace(0, 5, 5)

    y1 = 1.5 - 2 * x
    y2 = 2 - 2 * x
    y3 = x * 1
    y4 = 1.5 * x

    ax.fill_between(x, y1, y2, alpha = 0.5,)
    ax.fill_between(x, y3, y4, alpha = 0.5,)

    fig.set_figwidth(12)    #  ширина и
    fig.set_figheight(6)    #  высота "Figure"
    fig.set_facecolor('floralwhite')
    ax.set_facecolor('seashell')

    for elem in iterations:
        if elem[0][0].is_EmptySet is None and elem[1][0].is_EmptySet is None:
            ax.add_patch(
                 pth.Rectangle((elem[0][0].inf, elem[1][0].inf), (elem[0][0].sup - elem[0][0].inf), (elem[1][0].sup - elem[1][0].inf),
                               linewidth=1, edgecolor='r', facecolor='none'))

    plt.xlabel("x1")
    plt.ylabel("x2")
    plt.show()

def plotRadius(iterations):
    radius = []
    iter = []
    for i in range(len(iterations)):
        radius.append(math.sqrt(((iterations[i][0][0].sup - iterations[i][0][0].inf) / 2) ** 2 + ((iterations[i][1][0].sup - iterations[i][1][0].inf) / 2) ** 2))
        iter.append(i)
    plt.xlabel("Номер итерации")
    plt.ylabel("Радиус")
    plt.semilogy(iter, radius)
    plt.show()

def plotDistance(iterations):
    distance = []
    iter = []
    end_x = (iterations[len(iterations) - 1][0][0].sup + iterations[len(iterations) - 1][0][0].inf) / 2
    end_y = (iterations[len(iterations) - 1][1][0].sup + iterations[len(iterations) - 1][1][0].inf) / 2
    for i in range(len(iterations)):
        x = (iterations[i][0][0].sup + iterations[i][0][0].inf) / 2 - end_x
        y = (iterations[i][1][0].sup + iterations[i][1][0].inf) / 2 - end_y
        distance.append(math.sqrt(x ** 2 + y ** 2))
        iter.append(i)
    plt.xlabel("Номер итерации")
    plt.ylabel("Расстояние")
    plt.semilogy(iter, distance)
    plt.show()
