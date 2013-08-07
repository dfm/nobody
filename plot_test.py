#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)
import matplotlib.pyplot as pl
import numpy as np

data = np.loadtxt("orbit.out")
pl.plot(data[:, 1], data[:, 2])
pl.plot(data[:, 7], data[:, 8])
pl.savefig("orbit.png")
