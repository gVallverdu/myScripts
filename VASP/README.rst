VASP scripts
============

Script concerning post or pre treatments with VASP.

makeKpoins.py
    Create a KPOINTS file for a band structure calculation. This script use
    methods of pymatgen in order to compute and select high symetry lines in
    the first brillouin zone.

anaStruct.py
    Make histograms of distances and angles between atoms in a crystal
    structure read from a POSCAR/CONTCAR VASP file.

bands_Cu.py
    Build a band diagram with s, p, d contribution mapped on a RGB color scale.
    This script use pymatgen.
