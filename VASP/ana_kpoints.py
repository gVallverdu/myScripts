#!/usr/bin/env python

import os

# use pymatgen objects
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.vasp.outputs import Oszicar

# header
print("#  k    a       b       c     alpha   beta  gamma" +
      "     Econv           SinglePts        DE")

# OSZICAR files
oszicar = [f for f in os.listdir("./") if "OSZICAR_" in f]
oszicar.sort(key=lambda x: int(x.strip("OSZICAR_")))

# CONTCAR files
contcar = [f for f in os.listdir("./") if "CONTCAR_" in f]
contcar.sort(key=lambda x: int(x.strip("CONTCAR_")))

save = 0.0

# boucle sur les points k
for cont, osz in zip(contcar, oszicar):
    i = int(cont.strip("CONTCAR_"))

    # chargement de la structure
    p = Poscar.from_file(cont)
    a, b, c = p.structure.lattice.abc
    alpha, beta, gamma = p.structure.lattice.angles

    # lecture du fichier OSZICAR
    o = Oszicar(osz)
    elast = o.final_energy
    singlePoint = o.ionic_steps[0]["E0"]

    # affiche les resultats
    print("%4d %7.3f %7.3f %7.3f %6.1f %6.1f %6.1f %12.7f     %12.7f %12.3f" %
          (i, a, b, c, alpha, beta, gamma, elast, singlePoint,
           (singlePoint - save) * 1000.))

    # save energy to compute convergence
    save = singlePoint
