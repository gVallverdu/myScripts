#!/usr/bin/env python3
# -*- coding=utf-8 -*-

"""
Read the CRYSTAL output file after a SLABCUT and build a POSCAR for the surface model.
The slab model is put at the center of the box (-c option) or at the bottom of
the box (-b option). Use -h or --help to get help.
"""

__author__ = "Germain Vallverdu"
__email__ = "germain.vallverdu@univ-pau.fr"
__licence__ = "GPL"

import os
import argparse
from pymatgen.core.sites import Site
from pymatgen.core.structure import Lattice, Structure
from pymatgen.io.vasp.inputs import Poscar

def get_args():
    """ get CRYSTAL output file names and options """

    def exist(f):
        """ 'Type' for argparse - checks that file exists but does not open """
        if not os.path.isfile(f):
            raise argparse.ArgumentTypeError("{0} does not exist".format(f))
        return f

    parser = argparse.ArgumentParser(description=__doc__)

    # CRYSTAL output file name is mandatory
    parser.add_argument("cryout",
                        help="path to the CRYSTAL output file containing the SLABCUT run",
                        nargs=1,
                        metavar="FILE.out",
                        type=exist)
    # options
    parser.add_argument("-o", "--poscar",
                        help="Name of POSCAR output file",
                        metavar="POSCAR",
                        default="POSCAR_slabcut.vasp")
    parser.add_argument("-z", "--zvide",
                        help="Vacuum height in angstrom",
                        metavar="ZVIDE",
                        default=12.,
                        type=float)
    parser.add_argument("-c", '--center', dest='center',
                        help="Put the slab at the center of the box.",
                        action='store_true')
    parser.add_argument("-b", '--bottom', dest='bottom',
                        help="Put the slab at the bottom of the box.",
                        action='store_true')

    return parser.parse_args()


def read_struct(filename):
    """
    Read lattice and atom positions from the crystal output file after a
    SLABCUT run.
    """
    with open(filename, "r") as f:
        lignes = f.readlines()

    sectionCoord = False
    sites = list()
    for i, ligne in enumerate(lignes):

        if "DEFINITION OF THE NEW LATTICE VECTORS" in ligne:
            while "LATTICE PARAMETERS  (ANGSTROM" not in ligne:
                i += 1
                ligne = lignes[i]
            params = [float(val) for val in lignes[i+2].split()]
            slab_lattice = Lattice.from_parameters(*params[0:3], *params[3:])

        if "COORDINATES OF THE ATOMS BELONGING TO THE SLAB" in ligne:
            sectionCoord = True
            continue

        if sectionCoord:
            # on saute la premiere ligne
            if "LAB" in ligne:
                continue

            # si ligne blanche on sort
            if ligne.strip() == "" :
                break

            z = int(ligne.split()[1])
            valeurs = [float(val) for val in ligne[56:].split()]
            sites.append(Site(z, valeurs))

    return slab_lattice, sites

if __name__ == "__main__":

    # read args from command line
    args = get_args()
    slab_lattice, sites = read_struct(args.cryout[0])

    # --------------------
    # zmin et zmax
    # --------------------
    zmin = min([site.z for site in sites])
    zmax = max([site.z for site in sites])
    print("zmin  = {:12.7f}".format(zmin))
    print("zmax  = {:12.7f}".format(zmax))
    print("dz    = {:12.7f}".format(zmax - zmin))

    # --------------------------------------------------------------
    # new value of the lattice parameter perpendcular to the surface
    # --------------------------------------------------------------
    newc = zmax - zmin + args.zvide
    print("new c = {:12.7f}\n".format(newc))

    a, b, c = slab_lattice.abc
    newLattice = Lattice.from_parameters(a=a, b=b, c=newc,
                                         alpha=90., beta=90., gamma=slab_lattice.gamma)

    # --------------------------------------------
    # decalage des coordonnees entre 0 et 1 dans (xy)
    # conversion de z en coordonnees reduites
    # --------------------------------------------
    print(args)
    newSites = list()
    for site in sites:
        coords = site.coords
        coords[0] += .5
        coords[1] += .5
        if args.center:
            coords[2] = coords[2] / newLattice.c + .5
        elif args.bottom:
            coords[2] = (coords[2] - zmin) / newLattice.c
        else:
            coords[2] = coords[2] / newLattice.c

        newSites.append(Site(site.specie, coords))

    # -----------------------------------
    # make a structure and print POSCAR
    # -----------------------------------
    struct = Structure(newLattice, [site.specie for site in newSites], [site.coords for site in newSites])
    struct.sort()
    Poscar(struct).write_file(args.poscar)
    print(struct)
