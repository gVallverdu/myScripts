VASP scripts
============

Script concerning post or pre treatments with VASP.

v
    Do simple operations to restart calculations with VASP and output some results. Some
    methods need pymatgen (`http://pymatgen.org`_) to be used. Execute v -h to see the 
    complete documentation.

scofield.py
    Python module in order to interpret valence band using DOS and cross
    section from Scofield paper.

getCharges
    Compute atomic charges from a Bader caclculations done with the bader
    program of the University of Texas at Austin.

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

encut.job
    Bash script to compute gamma point energie as a function of the plane waves cutoff.

chgsum.f90
    Fortran program in order to sum up to CHGCAR files. This fortran program is
    two times faster than a python script.
