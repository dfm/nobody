#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, print_function, absolute_import

__all__ = ["NBodySystem"]

import numpy as np
from . import extrapolate


class NBodySystem(object):

    def __init__(self):
        self.masses = np.array([])
        self.positions = np.empty([0, 3])
        self.velocities = np.empty([0, 3])

    def add_particle(self, mass, position, velocity):
        self.masses = np.append(self.masses, mass)
        self.positions = np.vstack([self.positions, position])
        self.velocities = np.vstack([self.velocities, velocity])

    @property
    def accelerations(self):
        acc = np.zeros_like(self.positions)
        for i in range(len(self.positions)):
            for j in range(i + 1, len(self.positions)):
                delta = self.positions[i] - self.positions[j]
                d2 = np.dot(delta, delta)
                delta /= d2 * np.sqrt(d2)
                acc[i] -= self.masses[j] * delta
                acc[j] += self.masses[i] * delta
        return acc

    @property
    def vector(self):
        return np.concatenate([self.positions, self.velocities], axis=1)

    @vector.setter
    def vector(self, value):
        self.positions = value[:, :3]
        self.velocities = value[:, 3:]

    @property
    def gradient(self):
        return np.concatenate([self.velocities, self.accelerations], axis=1)

    @property
    def energy(self):
        d = np.eye(len(self.positions), dtype=bool)
        delta = self.positions[:, None, :] - self.positions[None, :, :]
        delta = np.sum(delta ** 2, axis=2)
        delta[d] = 1.0
        delta = self.masses[:, None] * self.masses[None, :] / np.sqrt(delta)
        delta[d] = 0.0
        ep = np.sum(delta[np.tri(len(self.positions), len(self.positions), -1,
                                 dtype=bool)])
        ek = np.sum(self.masses * np.sum(self.velocities ** 2, axis=1))
        return 0.5 * ek - ep

    def euler_step(self, h):
        self.vector = self.vector + h * self.gradient

    def leapfrog_step(self, h):
        self.velocities += h * self.accelerations
        self.positions += h * self.velocities

    def rk4_step(self, h):
        v0 = self.vector
        k1 = self.gradient
        self.vector = v0 + 0.5 * h * k1
        k2 = self.gradient
        self.vector = v0 + 0.5 * h * k2
        k3 = self.gradient
        self.vector = v0 + h * k3
        k4 = self.gradient
        self.vector = v0 + (k1+2*(k2+k3)+k4) * h / 6.0

    def midpoint_step(self, H, N=2):
        h = H / N
        x0 = self.vector
        x1 = x0 + h * self.gradient
        for n in range(2, N + 1):
            self.vector = x1
            x0, x1 = x1, x0 + 2 * h * self.gradient
        self.vector = x1
        self.vector = 0.5 * (x0 + x1 + h * self.gradient)

    def bs_step(self, H, eps=1e-8):
        initial = np.array(self.vector)
        ks = 2 * np.arange(1, 8 + 1)
        x = (H / ks) ** 2
        y = np.empty([len(x)] + list(initial.shape))
        tableau = np.empty_like(y)

        for i, k in enumerate(ks):
            self.vector = initial
            self.midpoint_step(H, N=k)
            y[i] = self.vector
            yout, yerr = extrapolate.poly_step(i, x, y, tableau)
            emax = np.max(np.abs(yerr))
            if emax < eps:
                break

        if emax > eps:
            self.vector = initial
            return self.bs_step(0.5 * H, eps=eps)

        self.vector = yout
        return H
