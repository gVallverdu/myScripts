#!/usr/bin/env python

import os

# use pymatgen objects
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.vasp.outputs import Oszicar

# header
print("# kspacing    grid     a       b       c       alpha beta    gamma"
      "  Econv           SinglePts        DE")


# OSZICAR files
oszicar = [f for f in os.listdir("./") if "OSZICAR_" in f]
oszicar.sort(key=lambda x: int(x.strip("OSZICAR_")))

# CONTCAR files
contcar = [f for f in os.listdir("./") if "CONTCAR_" in f]
contcar.sort(key=lambda x: int(x.strip("CONTCAR_")))

# OUTCAR files
outcar = [f for f in os.listdir("./") if "OUTCAR_" in f]
outcar.sort(key=lambda x: int(x.strip("OUTCAR_")))

kspacing = [float(val) for val in "1.4 0.9 0.8 0.6 0.45 0.4 0.35 0.3 0.275 0.25 0.225 0.2 0.19 0.17".split()]

save = 0.0

# boucle sur les points k
for cont, osz, out, ksp in zip(contcar, oszicar, outcar, kspacing):
    i = int(cont.strip("CONTCAR_"))

    # chargement de la structure
    p = Poscar.from_file(cont)
    a, b, c = p.structure.lattice.abc
    alpha, beta, gamma = p.structure.lattice.angles

    # lecture du fichier OSZICAR
    o = Oszicar(osz)
    elast = o.final_energy
    singlePoint = o.ionic_steps[0]["E0"]

    # lecture de l'OUTCAR
    with open(out, "r") as f:
        for line in f:
            if "generate k-points for:" in line:
                kx, ky, kz = [int(k) for k in line.split()[-3:]]
                break

    # affiche les resultats
    print("%10.2f  %2d %2d %2d %7.3f %7.3f %7.3f %6.1f %6.1f %6.1f %12.7f     %12.7f %12.3f" %
          (ksp, kx, ky, kz, a, b, c, alpha, beta, gamma, elast, singlePoint,
           (singlePoint - save) * 1000.))

    # save energy to compute convergence
    save = singlePoint
