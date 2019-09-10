#!/usr/bin/env python3

import os
import numpy as np
from pymatgen.io.vasp.outputs import Outcar
import matplotlib.pyplot as plt

ao = "1s"
zvals = [0.0, 0.2, 0.4, 0.5, 0.6, 0.8, 1.0]

foutcars = ["CL_%3.1f/OUTCAR" % z for z in zvals]
print(foutcars)

clz = list()
clz_e = list()

lines = "# Clz eps_i(1s) for all atoms\n"
for zval, foutcar in zip(zvals, foutcars):
    # print("OUTCAR = ", foutcar)
    out = Outcar(foutcar)
    cl = out.read_core_state_eigen()
    clz.append(cl[0][ao][-1])
    clz_e.append(cl[0][ao][-1] + out.efermi)
    line = "%5f" % zval
    for cl_at in cl:
        line += "%12.4f" % cl_at[ao][-1]
    print("%4.1f %10.4f %10.4f %10.4f" % (zval, out.efermi, cl[0][ao][-1], cl[0][ao][-1] + out.efermi))

    lines += line + "\n"
    #print(line)

#print(clz)
with open("cl.dat", "w") as f:
    f.write(lines)

params = np.polyfit(x=zvals, y=clz, deg=2)
a, b, c = params
print(params)
params = np.polyfit(x=zvals, y=clz_e, deg=2)
a, b, c = params
print(params)
fit = lambda x: a * x**2 + b * x + c

plt.plot(zvals, clz, "ro--", label="eps")
plt.plot(zvals, clz_e, "yo--", label="epsi - eF")
x = np.linspace(0, 1, 200)
y = fit(x)
plt.plot(x, y, "b-")
plt.plot((x[0], x[-1]), (y[0], y[-1]), "k-")
plt.legend()
plt.show()

print("integration:" ,np.trapz(y, x))
for z, e in zip(zvals, clz):
    print("%4.1f %10.4f" % (z, e))

