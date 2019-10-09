#!/usr/bin/env python

import numpy as np

# use pymatgen objects
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.vasp.outputs import Oszicar, Outcar

s = ["%.2f" % si for si in [0.9, 0.92, 0.94, 0.96, 0.98, 1.05]]
U = ["%.1f" % ui for ui in np.arange(1, 5.2, .2)]

lines = "# U  S  <mu_B>  s_mu_B  total mag  <d>  s_d  distances\n"
for ui in U:
    for si in s:
        folder = f"U{ui}/{si}"
        print(folder)

        out = Outcar(folder + "/OUTCAR")

        p = Poscar.from_file(folder + "/CONTCAR")
        structure = p.structure

        distances = list()
        mag_moms = list()
        for site, mag in zip(structure, out.magnetization):
            if site.symbol == "Co":
                mag_moms.append(mag)

                neighbors = structure.get_neighbors(site, radius)
                distance = [d for (_, d) in neighbors]
                distances.append(distance)

        distances = np.array(distances)
        ave_dist = distances.mean()
        std_dist = distances.sdt()

        mag_moms = np.array(mag_moms)
        ave_mag = mag_moms.mean()
        std_mag = mag_moms.std()

        osz = Oszicar(folder + "/OSZICAR")
        total_mag = osz.ionic_steps[-1]["mag"]


        lines += f"{ui:8.2f} {si:8.2f} "
        lines += f"{ave_mag:10.4f} {std_mag:10.4f} "
        lines += f"{total_mag:10.4f} "
        lines += f"{ave_dist:10.4f} {std_dist:10.4f} "
        lines += [f"{d:6.3f}" for d in distances]
        lines += "\n"
