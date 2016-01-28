#!/usr/bin/env python

import matplotlib.pyplot as pp
import numpy as np

size = (5, 1)
decomps = [(1, 1), (5, 1), (15, 1), (15, 2), (15, 4)]

pp.figure(figsize=size, dpi=80)
pp.axis("off")
pp.xlim(0, size[0])
pp.ylim(0, size[1])

for d in decomps:
    X = np.linspace(0, size[0], d[0] + 1)
    Y = np.linspace(0, size[1], d[1] + 1)

    for i in range(len(X) - 1):
        for j in range(len(Y) - 1):
            pp.axhspan(Y[j],
                       Y[j + 1],
                       X[i] / size[0],
                       X[i + 1] / size[1],
                       facecolor="teal",
                       edgecolor="1",
                       linewidth=2)

    pp.savefig("{}x{}.png".format(d[0], d[1]))
