#!/usr/bin/env python
# -*- coding=utf-8 -*-

"""
Create a KPOINTS file for a band structure calculation. This script use methods
of pymatgen in order to compute and select high symetry lines in the first
brillouin zone.
"""

__author__ = "Germain Salvato-Vallverdu"
__email__ = "germain.vallverdu@univ-pau.fr"
__licence__ = "GPL"
__date__ = "April 2014"

import pymatgen as mg
from pymatgen.symmetry.finder import SymmetryFinder
from pymatgen.symmetry.bandstructure import HighSymmKpath
from pymatgen.io.vaspio.vasp_input import Kpoints

# read structure
struct = mg.read_structure("POSCAR_NaPO3_A_uniq_a.vasp")

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
ndiv = 20
Kpoints.automatic_linemode(ndiv, ibz).write_file("KPOINTS")
