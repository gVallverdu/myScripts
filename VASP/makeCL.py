#!/usr/bin/env python
# -*- coding=utf-8 -*-

import os
import shutil
import pymatgen as mg
from pymatgen.io.vaspio.vasp_input import Poscar, Potcar


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

    line = "\nCore Level\n"
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
    poscar = Poscar.from_file(fposcar, read_velocities=False)

    # move up the selected atom and ad "_" to differentiate it
    site = poscar.structure.pop(iat)
    poscar.structure.insert(0, mg.DummySpecie("_" + site.specie.symbol), site.frac_coords)

    poscar.write_file("POSCAR_CL.vasp")

    # remove the _ symbol on 6th line
    with open("POSCAR_CL.vasp", "r") as f:
        lines = f.readlines()
        lines[5] = "".join(["%4s" % el.strip("_") for el in lines[5].split()]) + "\n"
        lines = "".join(lines)
    with open("POSCAR_CL.vasp", "w") as f:
        f.write(lines)


def makeRun(ao, job, clz=.5, poscar="POSCAR_CL.vasp", potcar="POTCAR_CL",
            incar="INCAR_CL", dirname=None):
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
        dirname (str): name of the job directory, if None "CL_%3.1f" % clz is used

    Return:
        dirname (str): name of the job directory
        jobname (str): name of the script for submission
    """

    # make job directory
    if not dirname:
        dirname = "CL_%3.1f" % clz
    dirname = os.path.join(os.getcwd(), dirname)
    os.mkdir(dirname)

    # copy files
    shutil.copy(poscar, os.path.join(dirname, "POSCAR"))
    shutil.copy(potcar, os.path.join(dirname, "POTCAR"))
    shutil.copy("KPOINTS", dirname)
    jobname = job + "_%3.1f" % clz
    shutil.copy(job, os.path.join(dirname, jobname))

    # complete and copy INCAR file
    with open(incar, "r") as fincar:
        incar_lines = fincar.read() + add_CL_tags("1s", clz)
    with open(os.path.join(dirname, "INCAR"), "w") as fincar:
        fincar.write(incar_lines)

    return dirname, jobname

if __name__ == "__main__":

    # pic de coeur Li1s
    buildPOTCAR("Li", "../opt/POTCAR")

    # make super cell
    struct = mg.Structure.from_file("../opt/CONTCAR")
    struct.make_supercell(4)

    p = mg.io.vasp.inputs.Poscar(struct)
    p.write_file("POSCAR_4.vasp")

    buildPOSCAR(1, "POSCAR_4.vasp")

    # make runs
    zval = [0., 0.2, 0.4, 0.5, 0.6, 0.8, 1.]
    print(zval)

    lines = "#!/bin/bash\n"
    for clz in zval:
        dirname, job = makeRun("1s", "jLi_bcc", clz)
        dirname = "/" + dirname.strip("/Users/gvallver/pyrene")
        print("cd " + dirname)
        print("sbatch " + job)
        lines += "cd %s \n" % dirname
        lines += "sbatch %s \n" % job

    with open("lancer.sh", "w") as f:
        f.write(lines)
