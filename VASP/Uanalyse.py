#!/usr/bin/env python

import numpy as np

# use pymatgen objects
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.vasp.outputs import Oszicar, Outcar

s = ["%.2f" % si for si in [0.9, 0.92, 0.94, 0.96, 0.98, 1.05]]
s = ["0.9", "0.92", "0.94", "0.96", "0.98", "1.05"]
U = ["%.1f" % ui for ui in np.arange(1, 5.2, .2)]

radius = 2.5

lines = "# U  S  <mu_B>  s_mu_B  total mag  <d>  s_d  distances\n"
for ui in U:
    for si in s:
        folder = f"U{ui}/{si}"
        print(folder)

        out = Outcar(folder + "/OUTCAR")
        magnetization = [data["tot"] for data in out.magnetization]
        

        p = Poscar.from_file(folder + "/CONTCAR")
        structure = p.structure

        distances = list()
        mag_moms = list()
        for site, mag in zip(structure, magnetization):
            if site.specie.symbol == "Co":
                mag_moms.append(mag)

                neighbors = structure.get_neighbors(site, radius)
                distances += [d for (_, d) in neighbors]

        distances = np.array(distances)
        ave_dist = distances.mean()
        std_dist = distances.std()

        mag_moms = np.array(mag_moms)
        ave_mag = mag_moms.mean()
        std_mag = mag_moms.std()

        osz = Oszicar(folder + "/OSZICAR")
        total_mag = osz.ionic_steps[-1]["mag"]


        lines += f"{ui:>8s} {si:>8s} "
        lines += f"{ave_mag:10.4f} {std_mag:10.4f} "
        lines += f"{total_mag:10.4f} "
        lines += f"{ave_dist:10.4f} {std_dist:10.4f} "
        lines += "".join([f"{d:6.3f}" for d in distances])
        lines += "\n"


with open("data.dat", "w") as f:
    f.write(lines)
