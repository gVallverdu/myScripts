#!/usr/bin/env python
# -*- coding=utf-8 -*-

"""
Create a KPOINTS file for a band structure calculation. This script use
methods of pymatgen in order to compute and select high symetry lines
in the first brillouin zone.

SYNTAX 
        makeKpoints.py [OPTIONS] [STRUCTURE FILE]

STRUCTURE FILE
        must contain a structure. For example a POSCAR file.

OPTIONS
        -d ndiv
                ndiv is a integer corresponding to the number of 
                k-points needed along each symetry line
"""

__author__ = "Germain Salvato-Vallverdu"
__email__ = "germain.vallverdu@univ-pau.fr"
__licence__ = "GPL"
__date__ = "April 2014"

import sys
import os

import pymatgen as mg
from pymatgen.symmetry.finder import SymmetryFinder
from pymatgen.symmetry.bandstructure import HighSymmKpath
from pymatgen.io.vaspio.vasp_input import Kpoints

# default args
fstruct = "POSCAR"
ndiv = 20

# read args
if len(sys.argv) > 1:
    if len(sys.argv) == 2:
        fstruct = sys.argv[1]
    else:
        fstruct = sys.argv[-1]
    if "-d" in sys.argv:
        print(sys.argv)
        try:
            ndiv = int(sys.argv[sys.argv.index("-d") + 1])
        except ValueError, IndexError:
            print("-d must be followed by an integer")
            exit(1)

# read structure
if os.path.exists(fstruct):
    struct = mg.read_structure(fstruct)
else:
    print("File %s does not exist" % fstruct)
    exit(1)

# symmetry information
struct_sym = SymmetryFinder(struct)
print("lattice type : {0}".format(struct_sym.get_lattice_type()))
print("space group  : {0} ({1})".format(struct_sym.get_spacegroup_symbol(),
                                     struct_sym.get_spacegroup_number()))

# Compute first brillouin zone
ibz = HighSymmKpath(struct)
ibz.get_kpath_plot(savefig="path.png")
print("ibz type     : {0}".format(ibz.name))

# print specific kpoints in the first brillouin zone
for key, val in ibz.kpath["kpoints"].items():
    print("%8s %s" % (key, str(val)))
 
# suggested path for the band structure
print("paths in first brillouin zone :")
for path in ibz.kpath["path"]:
    print(path)

# write the KPOINTS file
Kpoints.automatic_linemode(ndiv, ibz).write_file("KPOINTS")
