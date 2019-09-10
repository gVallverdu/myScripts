#!/usr/bin/env python3
# -*- coding=utf-8 -*-

"""
cry2cif\n\n

Read the last geometry corresponding to the CRYSTALLOGRAPHIC CELL on a
CRYSTAL09 output file and print it in a cif format. If geometry
optimization did not converge, input geometry is printed instead.
"""

# TODO:
#  * returns coordinates instead of write a file
#  * make functions for various formats

__author__ = "Germain Vallverdu"
__email__ = "germain.vallverdu@univ-pau.fr"
__licence__ = "GPL"

import os
import argparse
from pymatgen import Structure, Lattice
from crystalio import CrystalOutfile


def get_options():
    """ get options from command lines """

    parser = argparse.ArgumentParser(prog="cry2cif", description=__doc__)

    # mandatory argument is CRYSTAL output filename
    parser.add_argument("filename",
                        help="CRYSTAL output file",
                        metavar="FILENAME",
                        type=str)

    # choose either cif or POSCAR format
    parser.add_argument("-t", "--to",
                        help="output format: either cif or VASP (POSCAR)",
                        metavar="format",
                        default="cif",
                        choices=("cif", "vasp"),
                        type=str)

    # center slab or nanotubes
    parser.add_argument("-i", "--center",
                        help="move the slab or nanotubes in the center of the box",
                        action="store_true",
                        dest="center",
                        default=False)

    parser.add_argument("-n", "--num_structure",
                        help="Structure number to be extracted (default, the last)",
                        metavar="N",
                        default=-1,
                        type=int)

    # sort atom along z or x for slab or nanotubes
    parser.add_argument("-z", "--sortz",
                        help="Sort atoms along z axis (for slabs)",
                        dest="sortz",
                        action="store_true",
                        default=False)
    parser.add_argument("-x", "--sortx",
                        help="Sort atoms along x axis (for nanotubes)",
                        dest="sortx",
                        action="store_true",
                        default=False)

    # in the case of slabs or nanotubes, you have to give a value for b or c
    parser.add_argument("-b",
                        help="lattice parameter b",
                        metavar="b",
                        default=50,
                        type=float)
    parser.add_argument("-c",
                        help="lattice parameter c",
                        metavar="c",
                        default=50,
                        type=float)

    return parser.parse_args()


def cry2cif(filename, to="cif", center=False, sortx=False, sortz=False,
            b_dum=50, c_dum=50, istruct=-1):
    """
    Read a CRYSTAL output file and return the structure in a cif or POSCAR format.

    Args:
        filename (str): crystal output filename
        to (str): 'cif' or 'vasp', format of the output file (default is cif)
        center (bool): if True, the slab or nanotube is translated to the center of
                       the box (default is False)
        sortx (bool): Nanotube : if True, atoms are sorted along x axes (default is False).
        sortz (bool): slab : if True, atoms are sorted along z axes (default is False).
        b_dum (float): dummy lattice paramters b in angstrom for nanotubes (default 50 A)
        c_dum (float): dummy lattice paramters c in angstrom for nanotubes and slabs (default 50 A)
        istruct (int): structure to be extracted

    """
    cryout = CrystalOutfile(filename)

    print("title      : ", cryout.title)
    if cryout.group:
        print("group      : ", cryout.group)

    # print("Number of structure read: ", len(cryout.structures))

    if istruct == -1:
        print("structure  : Final structure")
        structure = cryout.final_structure
    else:
        print("structure  : Structure %d" % istruct)
        structure = cryout.get_structure(istruct)
    print("# atom     : ", structure.num_sites)
    print("composition: ", structure.composition)
    print("Cell parameters:")
    print("a     : %10.4f" % structure.lattice.a)
    print("b     : %10.4f" % structure.lattice.b)
    print("c     : %10.4f" % structure.lattice.c)
    print("alpha : %10.4f" % structure.lattice.alpha)
    print("beta  : %10.4f" % structure.lattice.beta)
    print("gamma : %10.4f" % structure.lattice.gamma)

    # ----------------------------------------------------------
    #  New b and c axes
    # ----------------------------------------------------------
    if cryout.slab:
        frac_coords = structure.frac_coords
        frac_coords[:, 2] *= structure.lattice.c / c_dum
        matrix = structure.lattice.matrix.copy()
        matrix[2, 2] = c_dum
        structure = Structure(Lattice(matrix), structure.species, frac_coords)
    if cryout.nanotube:
        frac_coords = structure.frac_coords
        frac_coords[:, 1] *= structure.lattice.c / c_dum
        frac_coords[:, 2] *= structure.lattice.b / b_dum
        matrix = structure.lattice.matrix.copy()
        matrix[1, 1] = b_dum
        matrix[2, 2] = c_dum
        structure = Structure(Lattice(matrix), structure.species, frac_coords)

    # ----------------------------------------------------------
    #  move slab or nanotube to the center of the box
    # ----------------------------------------------------------
    if center:
        if cryout.slab:
            coords = structure.frac_coords.copy()
            coords[:, 2] += .5
            structure = Structure(structure.lattice, structure.species, coords)
        elif cryout.nanotube:
            coords = structure.frac_coords
            coords += .5
            structure = Structure(structure.lattice, structure.species, coords)

    # ----------------------------------------------------------
    #  sort atom along x or z axis for slab
    # ----------------------------------------------------------
    if sortz:
        isort = 2
    elif sortx:
        isort = 0

    axes = {2: "z", 0: "x"}
    if sortz or sortx:
        print("\nSort atoms along %s" % axes[isort])
        data = zip(structure.species, structure.frac_coords)
        data = sorted(data, key=lambda d: d[-1][isort], reverse=True)

        species = [d[0] for d in data]
        coords = [d[1] for d in data]
        structure = Structure(structure.lattice, species, coords)

    # ----------------------------------------------------------
    # export in the given format
    # ----------------------------------------------------------
    basename, _ = os.path.splitext(filename)
    if to.lower() == "cif":
        ext = ".cif"
    elif to.lower() == "vasp":
        to = "POSCAR"
        ext = ".vasp"
    else:
        to = "POSCAR"
        ext = ".vasp"
    structure.to(to, filename=basename + ext)


if __name__ == "__main__":
    # get arguments
    args = vars(get_options())
    # rename some args
    args["b_dum"] = args.pop("b")
    args["c_dum"] = args.pop("c")
    args["istruct"] = args.pop("num_structure")
    # call main program
    cry2cif(**args)
