import numpy as np
import matplotlib.pyplot as pl
from nobody.nobody import NBodySystem

s1 = NBodySystem()
s1.add_particle(1.0, [1.0, 0, 0], [0, 0.5, 0])
s1.add_particle(2.0, [0.0, 0, 0], [0, -0.1, 0])

s2 = NBodySystem()
s2.add_particle(1.0, [1.0, 0, 0], [0, 0.5, 0])
s2.add_particle(2.0, [0.0, 0, 0], [0, -0.1, 0])

dt = 0.001
factor = 20
N = 10000
pos = np.zeros((N, 2, 3))
pos2 = np.zeros((factor * N, 2, 3))
energy = np.empty(N)
for i in range(N):
    h = s1.midpoint_step(dt, N=12)
    pos[i] = s1.positions
    energy[i] = s1.energy
    # for j in range(factor):
    #     s2.rk4_step(h / factor)
    #     pos2[factor*i + j] = s2.positions
    # print(h, np.sum((s2.positions - pos[i]) ** 2))

pl.plot(pos[:, 0, 0], pos[:, 0, 1])
pl.plot(pos[:, 1, 0], pos[:, 1, 1])
# pl.plot(pos2[:, 0, 0], pos2[:, 0, 1])
# pl.plot(pos2[:, 1, 0], pos2[:, 1, 1])
pl.savefig("euler.png")

pl.figure()
pl.plot(energy)
pl.savefig("energy.png")
