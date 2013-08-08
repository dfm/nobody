#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)
import matplotlib.pyplot as pl
import numpy as np

data = np.loadtxt("orbit.out")
pl.figure(figsize=(8, 8))
pl.plot(data[:, 1], data[:, 2], "+")
pl.plot(data[:, 7], data[:, 8], ".")
pl.plot(data[:, 13], data[:, 14])
pl.xlim(-250, 250)
pl.ylim(-250, 250)
pl.savefig("orbit.png")

pl.clf()
pl.plot(data[:, 0], data[:, 4], "k")
pl.plot(data[:, 0], -data[:, 16] * 1e-5, "--r")
pl.savefig("rv.png")
