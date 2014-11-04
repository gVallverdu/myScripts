#!/usr/bin/env python
# -*- coding=utf-8 -*-

"""
Example of DOS modulation
Apply moduleDOS module
"""

from __future__ import division, print_function

import moduleDOS
from pymatgen.io.vaspio.vasp_output import Vasprun

__author__ = "Germain Salvato-Vallverdu"
__email__ = "germain.vallverdu@univ-pau.fr"
__licence__ = "GPL"
__date__ = "Nov. 2014"

# options:
xmlfile = "dos.xml"
tofile = False

# print cross sections:
print("Cross sections :")
run = Vasprun(xmlfile, parse_dos=False, parse_eigen=False)
for el in run.initial_structure.composition.elements:
    print(moduleDOS.CrossSec.from_string(el.symbol))

# molate DOS:
moduleDOS.modulate(xmlfile, tofile)



