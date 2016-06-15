#!/usr/bin/env python3
# -*- coding=utf-8 -*-

"""
cry2cif

NAME
        cry2cif - extract geometry from a CRYSTAL09 output file

SYNTAX
        cry2cif [OPTIONS] ... [FILE]

DESCRIPTION
        Read the last geometry corresponding to the CRYSTALLOGRAPHIC CELL on a
        CRYSTAL09 output file and print it in a cif format. If geometry
        optimization did not converge, input geometry is printed instead.

        -h, --help
            print this help and exit.

AUTHOR
        Germain Vallverdu <germain.vallverdu@univ-pau.fr>
"""

__author__ = "Germain Vallverdu"
__email__ = "germain.vallverdu@univ-pau.fr"
__licence__ = "GPL"

import sys
import os
import datetime

def usage(code):
    """ cry2cif usage """
    print(__doc__)
    exit(code)

def cry2cif(args):
    """ lit la geometrie sur le fichier out de crystal et ecrit un fichier cif. La
    geometrie utilisee est la crystallographic cell """

    # ----------------------------------------------------------
    # gestion des options
    # ----------------------------------------------------------
    narg = len(args)
    if narg == 2:
        if args[1] == "-h" or args[1] == "--help":
            usage(0)
        else:
            if not os.path.exists(args[1]):
                print("Error : file {0} does not exist".format(args[1]))
                exit(1)
            else:
                outfilename = args[1].strip()
    else:
        print(args)
        print("Error : bad arguments")
        usage(1)


    # ----------------------------------------------------------
    # lecture du fichier output
    # ----------------------------------------------------------
    slab = False
    locGroupe = "SPACE GROUP"
    primitive = True
    locMaille = " PRIMITIVE CELL"
    with open(outfilename, "r") as f:
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
                print("nom      : {0}".format(phasename))

            elif locGroupe in line:
                group = line.split(":")[1].strip()
                print("groupe   : {0}".format(group))

            elif "TRANSFORMATION MATRIX PRIMITIVE-CRYSTALLOGRAPHIC CELL" in line:
                primitive = False
                locMaille = " CRYSTALLOGRAPHIC CELL "

            elif "FINAL OPTIMIZED GEOMETRY" in line:
                end = False
                break

        if end:
            print("WARNING:")
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

        if not slab:
            volume = float(line.split("=")[1].split()[0].strip(")"))
        f.readline()

        a, b, c, alpha, beta, gamma = f.readline().split()
        if slab:
            c = float(input("SLAB : c value ? : "))

        for i in range(4):
            f.readline()

        nom = list()
        red = list()
        Z = list()
        line = f.readline()
        while line != "\n":
            i, p, Zi, nomi, xi, yi, zi = line.split()
            if p == "F":
                # only one homologue is recorded
                line = f.readline()
                continue

            nom.append(nomi)
            Z.append(int(Zi))
            if slab:
                if float(zi) > c / 2.:
                    print("ERROR zi > c/2")
                    print("zi = " + zi)
                    print("c = " + str(c/2))
                    exit(1)
                red.append([float(xi), float(yi), float(zi) / float(c) ])
            else:
                red.append([float(xi), float(yi), float(zi)])

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
    composition = {nom[0]:1}
    for n in nom[1:]:
        if n in composition.keys():
            composition[n] += 1
        else:
            composition[n] = 1
    compo = ""
    for X, n in composition.items():
        compo += "{0}{1} ".format(X, n)

    print("compos   : {0}".format(compo.strip()))
    print("nat      : {0}".format(len(nom)))

    # ----------------------------------------------------------
    # write cif file
    # ----------------------------------------------------------
    cifname = os.path.splitext(outfilename)[0] + ".cif"
    lines =  "#---------------------------------------------------------------------\n"
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
    if not slab:
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
    for label, r in zip(nom, red):
        lines += "{:4s} {:20.12f} {:20.12f} {:20.12f} {:10.6f}\n".format(label, r[0], r[1], r[2], 1.)

    with open(cifname, "w") as f:
        f.write(lines)

if __name__ == "__main__":
    cry2cif(sys.argv)
