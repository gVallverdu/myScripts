#!/usr/bin/python
# -*- coding=utf-8 -*-

import sys
import pymatgen as mg
import numpy as np
from pymatgen.io.vaspio.vasp_input import Poscar

def look4VacantPosition():
    """ """

    s = mg.Structure.from_file(sys.argv[1])

    fill_s = s.copy()
    grid = mg.Structure(s.lattice, [], [])

    # grid step in frac coords
    step = 0.1
    # radius for looking for neighbors
    radius = 2.5
    # threshold to exclude vacant position
    els = s.composition.elements
    th = {els[0]: 1.7, els[1]: 2., els[2]: 2.}

    # grid ranges
    xmin = 0
    xmax = 1
    ymin = 0
    ymax = 1
    zmin = 0
    zmax = 1

    print("Parameters\n" + 10 * "-")
    print("Grid step          : %6.2f" % step)
    print("Neighbors radius   : %6.2f" % radius)
    print("threshold tolerence:")
    for el, dmax in th.items():
        print("X - %2s : %8.3f" % (el.symbol, dmax))

    vaclist = list()
    for xi in np.arange(xmin, xmax, step):
        for yi in np.arange(ymin, ymax, step):
            for zi in np.arange(zmin, zmax, step):
                vacsite = mg.PeriodicSite("X", [xi, yi, zi], s.lattice)
                neigh = s.get_neighbors(vacsite, radius)

                grid.append("X", [xi, yi, zi])

                add = True
                nLi = 0
                for neighsite, d in neigh:
                    dmax = th[neighsite.specie]
                    if d < dmax:
                        add = False
                        break
                    if neighsite.specie.symbol == "Li":
                        nLi += 1

                if nLi >= 2:
                    print("nLi >= 2")
                    add = False

                if add:
                    fill_s.append(vacsite.specie, vacsite.frac_coords)
                    vaclist.append([xi, yi, zi])

    print("\nResults\n" + 7 * "-")
    print("Write file 'POSCAR_fill.vasp': initial structure + vacant positions")
    Poscar(fill_s).write_file("POSCAR_fill.vasp")
    print("Write file 'POSCAR_grid.vasp': grid positions (%d atoms)" % len(grid))
    Poscar(grid).write_file("POSCAR_grid.vasp")
    print("\nNumber of positions : %d" % len(vaclist))
    if len(vaclist) < 30:
        for xyz in vaclist:
            print((3 * "%8.3f") % tuple(xyz))

if __name__ == "__main__":
    look4VacantPosition()

