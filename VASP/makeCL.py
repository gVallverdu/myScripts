#!/usr/bin/env python
# -*- coding=utf-8 -*-

from __future__ import division, print_function
import os
import shutil

def qnumber(oa):
    """ 
    return n and l quantum number from an atomic orbital defined as a string.

    args:
        oa (str): atomic orbital such as '1s', '3p' ...

    return:
        n, l: quantum number
    """
    oa = oa.strip()
    try:
        n = int(oa[0])
    except ValueError:
        print("Reading n: Error in AO name")
        print("OA = " + oa)
        exit(1)

    try:
        l = ["s", "p", "d", "f"].index(oa[1])
    except ValueError:
        print("Reading l: Error in AO name")
        print("OA = " + oa)
        exit(1)

    return n, l

def add_CL_tags(oa, clz=.5):
    """ 
    print core level part of INCAR file. The atom for which calculation is done, is
    supposed to be the first one (CLNT = 1).
    
    args:
        oa (str): atomic orbital such as '1s', '3p' ...
        clz (float): CLZ value
    """

    n, l = qnumber(oa)

    line = "Core Level\n"
    line += "  ICORELEVEL = 2\n"
    line += "  CLNT = 1\n"
    line += "  LVTOT = .TRUE.\n"
    line += "  CLN = %d\n" % n
    line += "  CLL = %d\n" % l
    line += "  CLZ = %3.1f\n\n" % clz

    return line

def buildPOTCAR(atname, fpotcar="POTCAR"):
    """
    Build POTCAR file according to atom for wich CL calculation is required. The POTCAR of
    the atom is placed in first postion.
    """
    from pymatgen.io.vaspio.vasp_input import Potcar
    potcar = Potcar.from_file(fpotcar)
    add = False
    for pot in potcar:
        if pot.symbol == atname:
            newpot = pot
            add = True
            break
    if not add:
        print("Atom " + atname + " not found in POTCAR")
        exit(1)
    else:
        potcar.insert(0, newpot)

    potcar.write_file("POTCAR_CL")

def buildPOSCAR(iat, fposcar="POSCAR"):
    """
    Build POSCAR file according to atom for wich CL calculation is required. The
    corresponding atom is placed in first position in the POSCAR.
    """
    from pymatgen.io.vaspio.vasp_input import Poscar
    poscar = Poscar.from_file(fposcar)

    site = poscar.structure.pop(iat)
    poscar.structure.insert(0, site.specie, site.frac_coords)
    poscar.write_file("POSCAR_CL.vasp")

def makeRun(ao, job, clz=.5, poscar="POSCAR_CL.vasp", potcar="POTCAR_CL", incar="INCAR_CL"):
    """
    Setup the job:

        1. Make job directory named CL_XX where XX is the clz value.
        2. Copy POSCAR, POTCAR and KPOINTS input files.
        3. Add Core Level tags to the INCAR file according to the AO.
        4. Copy job file in the job directory and add clz value to the name.


    Note that POTCAR and POSCAR file have to be previously checked. The atom for which the
    corelevel calculation will be done is the first one. See buildPOTCAR and buildPOSCAR
    functions.
    
    Args:
        oa (str): atomic orbital such as '1s', '3p' ...
        job (str): name of job file to submit calculation.
        clz (float): value of CLZ INCAR tag
        poscar (str): name of POSCAR file.
        potcar (str): name of POTCAR file.
        incar (str): name of the initial INCAR file.

    Return:
        dirname (str): name of the job directory
        jobname (str): name of the script for submission
    """

    # make job directory
    dirname = os.path.join(os.getcwd(), "CL_%3.1f" % clz)
    os.mkdir(dirname)

    # copy files
    shutil.copy(poscar, os.path.join(dirname, "POSCAR"))
    shutil.copy(potcar, os.path.join(dirname, "POTCAR"))
    shutil.copy("KPOINTS", dirname)
    jobname = job + "_%3.1f" % clz
    shutil.copy(job, os.path.join(dirname, jobname))

    # complete and copy INCAR file
    incar = open(incar, "r").read() + add_CL_tags("1s", clz)
    open(dirname + "/INCAR", "w").write(incar)

    return dirname, jobname

if __name__ == "__main__":

    # pic de coeur O1s
    buildPOTCAR("O")

    # O71 est le plus proche du Si central
    buildPOSCAR(71, "CONTCAR")

    zval = [x/10 for x in range(11)]
    print(zval)

    for clz in zval:
        dir, job = makeRun("1s", "jSiO2", clz)
        print("cd " + dir)
        print("qsub " + job)
