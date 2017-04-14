#!/usr/bin/env python
# -*- coding=utf-8 -*-

from __future__ import print_function, division

"""
script anaStruct
----------------

L'objectif de ce script est de calculer tous les angles et toutes les liaisons
à partir d'un fichier vasprun.xml ou POSCAR/CONTCAR.

"""

import os
import sys
import numpy as np
from scipy.stats import norm
import pymatgen as mg
import matplotlib.pyplot as plt


def anaStruct(poscar, calcDistance=True, calcAngle=True, plot=True):
    """ main program """

    # load data
    if not os.path.exists(poscar):
        print("file %s does not exist" % poscar)
        exit(1)
    else:
        print("Read file ", poscar)
        struct = mg.Structure.from_file(poscar)

    # distance analysis
    sigma = .005
    npts = 500
    xmin = .9
    xmax = 2.
    if calcDistance:
        x, data = getDistances(struct, sigma, npts, xmin, xmax)

        if plot:
            for bond, h in data.items():
                plt.plot(x, h, label=bond)
            plt.xlabel(r"distance   /   $\AA$")
            plt.ylabel("Histogram")
            plt.grid()
            plt.legend()
            plt.title("Histogram of distances")
            plt.show()

    # angle analysis
    amin = 80.
    amax = 100.
    cutoff = 3.0
    centralAtom = "Mn"
    ligandAtom = "O"
    sigma = .5
    npts = 200
    if calcAngle:
        x, h = getAngles(struct, sigma, npts, amin, amax, cutoff, centralAtom, ligandAtom)

        if plot:
            plt.plot(x, h, "k-")
            plt.xlabel("Angle   /   degree")
            plt.ylabel("Histogram")
            plt.grid()
            plt.title("Histogram of angles %s-%s-%s" % (ligandAtom,
                                                        centralAtom, ligandAtom))
            plt.show()


def getAngles(struct, sigma, npts, amin, amax, cutoff, centralAtom, ligandAtom):
    """ compute all angles  """

    # info
    print("\nBending angles analyses :")
    print("-------------------------")
    print("Central atom : %s" % centralAtom)
    print("Ligand atom  : %s" % ligandAtom)
    print("cutoff       : %f" % cutoff)

    # central sites
    central_sites = [site for site in struct
                     if site.specie == mg.Element(centralAtom)]

    # x values for histogram
    aval = np.linspace(amin, amax, npts)

    # data
    data = np.zeros(npts)

    print("\nCoordinence")
    for csite in central_sites:
        neighbors = [(site, d) for site, d in struct.get_neighbors(csite, cutoff)
                     if site.specie == mg.Element(ligandAtom)]
        iat = struct.index(csite)
        print("%5s(%2d)  : %d" % (csite.specie.symbol, iat, len(neighbors)))

        #neighbors.append((distance, xij))

        # compute bending angles and build histogram
        for i in range(len(neighbors)):
            isite, di = neighbors[i]
            xic = isite.coords - csite.coords
            for j in range(i + 1, len(neighbors)):
                jsite, dj = neighbors[j]
                xjc = jsite.coords - csite.coords
                scal = np.dot(xic, xjc) / di / dj
                if np.fabs(scal) > 1.:
                    print("Dot product error : %f" % scal)
                    raise ValueError
                angle = np.degrees(np.arccos(scal))

                # add a gaussian function
                data += norm.pdf(aval, loc=angle, scale=sigma)

    # print data
    line = "# structural analysis\n"
    line += "# column 1 : angle (degree)\n"
    line += "# column 2 : histogram of angle %s-%s-%s\n" % (ligandAtom, centralAtom, ligandAtom)
    for a, h in zip(aval, data):
        line += "%12.4f %12.4f\n" % (a, h)

    print("\n=> Print file anaAngles.dat\n")
    open("anaAngles.dat", "w").write(line)

    return aval, data


def getDistances(struct, sigma, npts, xmin, xmax):
    """ compute all distances """

    # system composition
    atomTypes = [e.symbol for e in struct.composition.elements]
    NAtomsType = len(struct.composition.elements)

    print("\nBond length analyses :")
    print("----------------------")
    print("Atoms   : ", " ; ".join(atomTypes))
    print("N atoms : ", NAtomsType)

    # bond types
    bondsList = list()
    for ityp in range(NAtomsType):
        for jtyp in range(ityp, NAtomsType):
            if atomTypes[ityp] > atomTypes[jtyp]:
                bondsList.append(atomTypes[ityp] + "-" + atomTypes[jtyp])
            else:
                bondsList.append(atomTypes[jtyp] + "-" + atomTypes[ityp])

    NBondsType = len(bondsList)
    print("N bonds : ", (NBondsType))
    print("Bonds   : " + " ; ".join(bondsList))

    # x values for histogram
    xval = np.linspace(xmin, xmax, npts)

    # histogram
    data = dict()

    # compute all distances
    Natoms = len(struct)
    for iat in range(Natoms):
        for jat in range(iat + 1, Natoms):

            # compute the cartesian coordinate and the distance
            distance = struct.get_distance(iat, jat)

            # xmax act as a cutoff
            if distance > xmax:
                continue

            # select the relevant bond
            if struct[iat].specie.symbol > struct[jat].specie.symbol:
                selectedBond = struct[iat].specie.symbol + "-" + struct[jat].specie.symbol
            else:
                selectedBond = struct[jat].specie.symbol + "-" + struct[iat].specie.symbol

            # add a gaussian
            if selectedBond in data:
                data[selectedBond] += norm.pdf(xval, loc=distance, scale=sigma)
            else:
                data[selectedBond] = norm.pdf(xval, loc=distance, scale=sigma)

    # print data
    line = "# structural analysis \n"
    line += "# column 1 : x (A)\n"
    for i, bond in enumerate(data.keys()):
        line += "# column %d : bond %s \n" % (i + 2, bond)
    for i, x in enumerate(xval):
        line += "%10.4f" % x
        for key in data.keys():
            line += " %10.4f" % data[key][i]
        line += "\n"

    print("\n=> Print file anaDistances.dat\n")
    open("anaDistances.dat", "w").write(line)

    return xval, data

if __name__ == "__main__":
    # poscar name
    if len(sys.argv) == 2:
        poscar = sys.argv[1]
    else:
        poscar = "CONTCAR"
    anaStruct(poscar, calcAngle=False, plot=True)
