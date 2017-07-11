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
import datetime
import argparse
from math import pi, cos, sin, sqrt


def usage(code):
    """ cry2cif usage """
    print(__doc__)
    exit(code)


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
                        choices=("cif", "POSCAR"),
                        type=str)

    # center slab or nanotubes
    parser.add_argument("-i", "--center",
                        help="move the slab or nanotubes in the center of the box",
                        action="store_true",
                        dest="center",
                        default=False)

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


def cry2cif(filename, filefmt="cif", center=False, sortx=False, sortz=False,
            b_dum=50, c_dum=50):
    """
    Read a CRYSTAL output file and return the structure in a cif or POSCAR format.

    Args:
        filename (str): crystal output filename
        filefmt (str): 'cif' or 'POSCAR', format of the output file (default is cif)
        center (bool): if True, the slab or nanotube is translated to the center of
                       the box (default is False)
        sortx (bool): Nanotube : if True, atoms are sorted along x axes (default is False).
        sortz (bool): slab : if True, atoms are sorted along z axes (default is False).
        b_dum (float): dummy lattice paramters b in angstrom for nanotubes (default 50 A)
        c_dum (float): dummy lattice paramters c in angstrom for nanotubes and slabs (default 50 A)

    """

    # ----------------------------------------------------------
    # lecture du fichier output
    # ----------------------------------------------------------
    slab = False
    nanotube = False
    locGroupe = "SPACE GROUP"
    primitive = True
    locMaille = " PRIMITIVE CELL"
    with open(filename, "r") as f:
        # read general data
        line = f.readline()
        end = True
        while line != "":
            line = f.readline()
            if "SLAB CALCULATION" in line:
                slab = True
                locGroupe = "PLANE GROUP"

            elif "SLABCUT" in line:
                slab = True

            elif "EEEEEEEEEE STARTING" in line:
                phasename = f.readline().strip()
                print("title    : {0}".format(phasename))

            elif locGroupe in line:
                group = line.split(":")[1].strip()
                print("group    : {0}".format(group))

            elif "TRANSFORMATION MATRIX PRIMITIVE-CRYSTALLOGRAPHIC CELL" in line:
                primitive = False
                locMaille = " CRYSTALLOGRAPHIC CELL "

            elif "FINAL OPTIMIZED GEOMETRY" in line:
                end = False
                break

            elif "SLAB GENERATED" in line:
                slab = True
                group = ""
                print("\nSLAB GENERATED")

            elif "CONSTRUCTION OF A NANOTUBE FROM A SLAB" in line:
                nanotube = True
                slab = False
                group = ""
                print("\nNANOTUBE FROM SLAB")

        if end:
            print("\nWARNING:")
            print("Optimisation did not converge, final optimized geometry not found.")
            print("Input geometry will be read instead.\n")
            f.seek(0)
            line = f.readline()
            while " GEOMETRY FOR WAVE FUNCTION " not in line:
                line = f.readline()

        # read geometry located at locMaille
        line = f.readline()
        while locMaille not in line:
            line = f.readline()

        if not slab and not nanotube:
            volume = float(line.split("=")[1].split()[0].strip(")"))
        f.readline()

        # lattice parameters
        a, b, c, alpha, beta, gamma = [float(val) for val in f.readline().split()]
        if slab:
            c = c_dum
        if nanotube:
            b = b_dum
            c = c_dum

        # if slab:
        #     c = float(input("\nSLAB : c value ? : "))
        #
        # if nanotube:
        #     b = float(input("\nNANOTUBE : b value ? : "))
        #     c = float(input("\nNANOTUBE : c value ? : "))

        for i in range(4):
            f.readline()

        # read coordinates
        nom = list()    # atom names
        red = list()    # reduce coordinates
        Z = list()      # atomic number
        uniq = list()   # True if atom belong to the asymmetric unit
        radius = list() # distance from the axes of the nanotube
        line = f.readline()
        while line != "\n":
            if nanotube:
                i, p, Zi, nomi, xi, yi, zi, ri = line.split()
            else:
                i, p, Zi, nomi, xi, yi, zi = line.split()

            xi = float(xi)
            yi = float(yi)
            zi = float(zi)

            Z.append(int(Zi))
            if nanotube:
                radius.append(float(ri))
                
            if p == "F":
                uniq.append(False)
            else:
                uniq.append(True)

            if len(nomi) == 2:
                nom.append(nomi[0] + nomi[1].lower())
            else:
                nom.append(nomi)

            if slab:
                if zi > c / 2.:
                    print("ERROR zi > c / 2")
                    print("zi = ", zi)
                    print("c = ", c / 2)
                    exit(1)
                red.append([xi, yi, zi / c])
            elif nanotube:
                if zi > c / 2.:
                    print("ERROR zi > c / 2")
                    print("zi = ", zi)
                    print("c = " , c / 2)
                    exit(1)
                if yi > b / 2.:
                    print("ERROR yi > b / 2")
                    print("yi = ", yi)
                    print("b = ", b / 2)
                    exit(1)
                red.append([xi, yi / b, zi / c])
            else:
                red.append([xi, yi, zi])

            line = f.readline()

    if end:
        print("cell     : guess geometry")
    else:
        if primitive:
            print("cell     : primitive")
        else:
            print("cell     : crystallographic")
    print("a        : {0}".format(a))
    print("b        : {0}".format(b))
    print("c        : {0}".format(c))
    print("alpha    : {0}".format(alpha))
    print("beta     : {0}".format(beta))
    print("gamma    : {0}".format(gamma))

    # ----------------------------------------------------------
    #  system composition
    # ----------------------------------------------------------
    composition = dict()
    atomTypes = list()
    for name in nom:
        if name in composition:
            composition[name] += 1
        else:
            composition[name] = 1
            atomTypes.append(name)
    compo = ""
    for X, n in composition.items():
        compo += "{0}_{1} ".format(X, n)

    print("compos   : {0}".format(compo.strip()))
    print("nat      : {0}".format(len(nom)))

    # ----------------------------------------------------------
    #  move slab or nanotube to the center of the box
    # ----------------------------------------------------------
    if center:
        if slab:
            red = [[r[0], r[1], r[2] + .5] for r in red]
        elif nanotube:
            red = [[r[0] + .5, r[1] + .5, r[2] + .5] for r in red]

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
        data = zip(nom, uniq, radius, red)
        dataSorted = sorted(data, key=lambda f:f[-1][isort], reverse=True)

        red_final = [ired for iname, iuniq, iradius, ired in dataSorted]
        uniq_final = [iuniq for iname, iuniq, iradius, ired in dataSorted]
        name_final = [iname for iname, iuniq, iradius, ired in dataSorted]
        radius_final = [iradius for iname, iuniq, iradius, ired in dataSorted]

    else:
        red_final = red
        uniq_final = uniq
        name_final = nom
        radius_final = radius

    if filefmt == "cif":
        # ----------------------------------------------------------
        # write cif file
        # ----------------------------------------------------------
        outname = os.path.splitext(filename)[0] + ".cif"
        lines = "#---------------------------------------------------------------------\n"
        lines += "# Date : {0}\n".format(datetime.datetime.now().strftime("%A %d %B %Y, %H:%M:%S"))
        lines += "# directory : {0}\n".format(os.getcwd())
        lines += "# hostname : {0}\n".format(os.uname()[1])
        lines += "#---------------------------------------------------------------------\n"
        lines += "#\n"
        lines += "# This file contains the last CRYSTALLOGRAPHIC CELL read on a CRYSTAL09\n"
        lines += "# output file and may be readable by a visualization tool such as VESTA :\n"
        lines += "# http://jp-minerals.org/vesta/en/\n"
        lines += "#\n"
        lines += "# Cell parameters\n"
        lines += "_pd_phase_name                    '{0}'\n".format(phasename)
        lines += "_cell_angle_alpha                 {0}\n".format(alpha)
        lines += "_cell_angle_beta                  {0}\n".format(beta)
        lines += "_cell_angle_gamma                 {0}\n".format(gamma)
        lines += "_cell_length_a                    {0}\n".format(a)
        lines += "_cell_length_b                    {0}\n".format(b)
        lines += "_cell_length_c                    {0}\n".format(c)
        if not slab and not nanotube:
            lines += "_cell_volume                      {0}\n".format(volume)
        lines += "_symmetry_space_group_name_H-M    '{0}'\n".format(group)
        lines += "_chemical_formula_sum             '{0}'\n".format(compo)

        lines += "\n"
        lines += "loop_\n"
        lines += "_atom_site_label\n"
        lines += "_atom_site_fract_x\n"
        lines += "_atom_site_fract_y\n"
        lines += "_atom_site_fract_z\n"
        lines += "_atom_site_occupancy\n"
        for name, p, r in zip(name_final, uniq_final, red_final):
            if p:
                lines += "%4s" % name
                lines += "%20.12f %20.12f %20.12f" % tuple(r)
                lines +=  "   1.\n"
            if not p and (slab or nanotube):
                lines += "%4s" % name
                lines += "%20.12f %20.12f %20.12f" % tuple(r)
                lines +=  "   1.\n"

        if nanotube:
            lines += "\n"
            lines += "# distance from the center of the nanotube (in A)\n"
            lines += "loop_\n"
            lines += "_atom_site_index\n"
            lines += "_atom_site_label\n"
            lines += "_atom_site_radius\n"
            for i, (name, radius) in enumerate(zip(name_final, radius_final)):
                lines += "%5d %4s %10.3f\n" % (i, name, radius)

    elif filefmt == "POSCAR":
        outname = "POSCAR_" + os.path.splitext(filename)[0] + ".vasp"
        lines = "Structure from %s\n" % filename
        lines += " 1.0\n"

        # compute lattice vectors
        alphar = alpha * pi / 180.0
        betar  = beta  * pi / 180.0
        gammar = gamma * pi / 180.0

        veca = [a, 0., 0.]
        vecb = [b * cos(gammar), b * sin(gammar), 0.]

        vecc = [0.,0.,0.]
        vecc[0] = c * cos( betar )
        cy = (cos(alphar) - cos(gammar) * cos(betar)) / sin(gammar)
        vecc[1] = c * cy
        cz = sqrt((sin(betar))**2 - cy**2)
        vecc[2] = c * cz

        lines += "%20.12f %20.12f %20.12f\n" % tuple(veca)
        lines += "%20.12f %20.12f %20.12f\n" % tuple(vecb)
        lines += "%20.12f %20.12f %20.12f\n" % tuple(vecc)

        # atom names and atom number
        lines += "".join(["%4s" % name for name in atomTypes]) + "\n"
        lines += "".join(["%4d" % composition[name] for name in atomTypes]) + "\n"

        # sort coordinates according to atom name
        redSorted = list()
        for name in atomTypes:
            for iname, r in zip(name_final, red_final):
                if iname == name:
                    redSorted.append(r)

        radius_out = "# RAW DATA FROM FILE %s\n" % filename
        radius_out += "# LATTICE PARAMETERS:\n"
        radius_out += "# a     = %20.12f\n" % a
        radius_out += "# b     = %20.12f\n" % b
        radius_out += "# c     = %20.12f\n" % c
        radius_out += "# alpha = %20.12f\n" % alpha
        radius_out += "# beta  = %20.12f\n" % beta
        radius_out += "# gamma = %20.12f\n" % gamma
        radius_out += "# LATTICE VECTORS:\n"
        radius_out += "# vector A = %20.12f %20.12f %20.12f\n" % tuple(veca)
        radius_out += "# vector B = %20.12f %20.12f %20.12f\n" % tuple(vecb)
        radius_out += "# vector C = %20.12f %20.12f %20.12f\n" % tuple(vecc)
        radius_out += "#\n# column 1: atom number\n"
        radius_out += "# column 2: atom name\n"
        radius_out += "# column 3: fractional coordinate x\n"
        radius_out += "# column 4: fractional coordinate y\n"
        radius_out += "# column 5: fractional coordinate z\n"
        radius_out += "# column 6: distance from the nanotube axes (angstrom)\n"
        for iat, (name, red, radius) in enumerate(zip(name_final, red_final, radius_final)):
            radius_out += "%4d %4s" % (iat, name)
            radius_out += "%20.12f %20.12f %20.12f" % tuple(red)
            radius_out += "%10.3f\n" % radius
        with open(os.path.splitext(filename)[0] + "_radius.dat", "w") as f:
            f.write(radius_out)

        # fractional coordinates
        lines += "Direct\n"
        for r in redSorted:
            lines += "%20.12f %20.12f %20.12f\n" % tuple(r)

    with open(outname, "w") as f:
        f.write(lines)

if __name__ == "__main__":
    # get arguments
    args = vars(get_options())
    filename = args["filename"]
    filefmt = args["to"]
    center = args["center"]
    sortx = args["sortx"]
    sortz = args["sortz"]
    b_dum = args["b"]
    c_dum = args["c"]

    # call main program
    cry2cif(filename, filefmt, center, sortx, sortz, b_dum, c_dum)
