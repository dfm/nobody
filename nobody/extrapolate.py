#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, print_function, absolute_import

__all__ = ["polynomial", "poly_step", "lagrange"]

import numpy as np


def polynomial(x, y):
    tableau = np.empty_like(y)
    for i in range(len(x)):
        yout, yerr = poly_step(i, x, y, tableau)
    return yout


def poly_step(i, x, y, tableau):
    yout = np.array(y[i])
    yerr = np.array(y[i])
    if i == 0:
        tableau[0] = np.array(y[0])
        return yout, yerr

    c = np.array(y[i])
    for k in range(1, i):
        delta = 1.0 / (x[i - k] - x[i])
        f1 = x[i] * delta
        f2 = x[i - k] * delta
        delta, tableau[k] = c - tableau[k], yerr
        yerr = f1 * delta
        c = f2 * delta
        yout += yerr
    tableau[i] = yerr
    return yout, yerr


def lagrange(x, y):
    result = 0.0
    for i, xi in enumerate(x):
        val = y[i]
        for j, xj in enumerate(x):
            if i != j:
                val *= xj / (xj - xi)
        result += val
    return result


if __name__ == "__main__":
    x = np.linspace(0, 1, 10)[:1:-1]
    y = np.cos(x)
    print(polynomial(x, y))
    print(lagrange(x, y))
