#!/usr/bin/python
# -*- coding=utf-8 -*-

"""
script anaStruct
----------------

L'objectif de ce script est de calculer tous les angles et toutes les liaisons
à partir d'un fichier vasprun.xml ou POSCAR/CONTCAR.

"""

import os
import numpy as np
import sys
from crystal import Crystal

# 27 cubes
trans = np.array([[0. , 0., 0.],
                  [ 1., 0., 0.],
                  [ 1., 1., 0.],
                  [ 1.,-1., 0.],
                  [ 1., 0., 1.],
                  [ 1., 0.,-1.],
                  [ 1., 1., 1.],
                  [ 1.,-1., 1.],
                  [ 1., 1.,-1.],
                  [ 1.,-1.,-1.],
                  [-1., 0., 0.],
                  [-1., 1., 0.],
                  [-1.,-1., 0.],
                  [-1., 0., 1.],
                  [-1., 0.,-1.],
                  [-1., 1., 1.],
                  [-1.,-1., 1.],
                  [-1., 1.,-1.], 
                  [-1.,-1.,-1.],
                  [ 0., 1., 0.],
                  [ 0.,-1., 0.],
                  [ 0., 0., 1.],
                  [ 0., 0.,-1.],
                  [ 0., 1., 1.],
                  [ 0.,-1., 1.],
                  [ 0., 1.,-1.],
                  [ 0.,-1.,-1.]])

def g(x, x0, sigma):
    """ return a gaussian function """
    return 1. / (sigma * np.sqrt(2. * np.pi)) * np.exp(-(x - x0)**2 / (2. * sigma**2))

# vectorize gaussian function
#vgauss = sp.vectorize(g, excluded = ["x0", "sigma"])

def anaStruct():
    """ main program """

    calcDistance = True
    calcAngle    = True

    # poscar name
    if len(sys.argv) == 2:
        poscar = sys.argv[1]
    else:
        poscar = "CONTCAR"

    # load data
    if not os.path.exists(poscar):
        print("file {0} does not exist".format(poscar))
        exit(1)
    else:
        print("Read file " + poscar)
        struct = Crystal.fromPOSCAR(poscar, verbose = False)

    # distance analysis
    sigma     = .01
    npts      = 300
    xmin      = 1.5
    xmax      = 2.5
    if calcDistance:
        getDistances(struct, sigma, npts, xmin, xmax)

    ## angle analysis
    amin        = 80.
    amax        = 100.
    cutoff      = 3.0
    centralAtom = "Co"
    ligandAtom  = "O"
    sigma       = .5
    npts        = 200
    if calcAngle:
        getAngles(struct, sigma, npts, amin, amax, cutoff, centralAtom, ligandAtom)

def getAngles(struct, sigma, npts, amin, amax, cutoff, centralAtom, ligandAtom):
    """ compute all angles  """

    # info
    print("\nBending angles analyses :")
    print(  "-------------------------")
    print("Central atom : %s" % centralAtom)
    print("Ligand atom  : %s" % ligandAtom)
    print("cutoff       : %f" % cutoff)

    # x values for histogram
    aval = np.linspace(amin, amax, npts)

    # data
    data = np.zeros(npts)
    nNeighbor = list()

    for iat in range(struct.Natoms):

        if struct.atomNames[iat] != centralAtom:
            continue

        riat = np.array([val for val in struct.redCoord[iat]]) # hard copy ?

        neighbors = list()
        for jat in range(struct.Natoms):

            if struct.atomNames[jat] != ligandAtom:
                continue

            # compute distance for each possible image
            rjat0 = np.array([val for val in struct.redCoord[jat]]) # hard copy ?
            for tr in trans:
                rjat = rjat0 + tr

                xij = struct.red2cart(rjat - riat)
                distance = np.sqrt(xij[0]**2 + xij[1]**2 + xij[2]**2)

                if distance > cutoff:
                    continue

                same = False
                for neigbhor in neighbors:
                    dn, xn = neigbhor
                    d = np.sqrt((xn[0] - xij[0])**2 \
                              + (xn[1] - xij[1])**2 \
                              + (xn[2] - xij[2])**2)
                    if d < 1.e-5:
                        same = True
                        break

                if not same:
                    neighbors.append((distance, xij))

        # compute bending angles and build histogram
        n = len(neighbors)
        nNeighbor.append((iat, struct.atomNames[iat], n))
        for i in range(n):
            di, xi = neighbors[i]
            for j in range(i + 1, n):
                dj, xj = neighbors[j]
 
                scal = np.dot(xi, xj) / di / dj
                if np.fabs(scal) > 1.:
                    print("Dot product error : %f" % scal)
                    exit(1)
                angle = np.arccos(scal) * 180.0 / np.pi

                # add a gaussian
                data += np.array([g(a, angle, sigma) for a in aval])

    print("\nCoordinence")
    for iat, name, n in nNeighbor:
        print("%5s(%2d)  : %d" % (name.rjust(5), iat, n))
    # print data
    line  = "# structural analysis %s\n" % struct.name
    line += "# column 1 : angle (degree)\n"
    line += "# column 2 : histogram of angle %s-%s-%s\n" % (ligandAtom, centralAtom, ligandAtom)
    for a, h in zip(aval, data):
        line += "%12.4f %12.4f\n" % (a, h)

    print("\n=> Print file anaAngles.dat\n")
    open("anaAngles.dat", "w").write(line)


def getDistances(struct, sigma, npts, xmin, xmax):
    """ compute all distances """

    # system composition
    atomTypes = list()
    for name in struct.atomNames:
        if name not in atomTypes:
            atomTypes.append(name.strip())

    NAtomsType = len(atomTypes)

    print("\nBond length analyses :")
    print(  "----------------------")
    print("Atoms   : {0}".format(" ; ".join(atomTypes)))
    print("N atoms : {0}".format(NAtomsType))

    # bond types
    bondsList = list()
    for ityp in range(NAtomsType):
        for jtyp in range(ityp + 1, NAtomsType):
            bondsList.append(atomTypes[ityp] + "-" + atomTypes[jtyp])

    NBondsType = len(bondsList)
    print("N bonds : {0}".format(NBondsType))
    print("Bonds   : {0}".format(" ; ".join(bondsList)))

    # x values for histogram
    xval = np.linspace(xmin, xmax, npts)

    # histogram
    data = dict()
    for bond in bondsList:
        data[bond] = np.zeros(npts)

    # compute all distances
    for iat in range(struct.Natoms):
        for jat in range(iat + 1, struct.Natoms):

            # continue if same atom type
            if struct.atomNames[iat] == struct.atomNames[jat]:
                continue

            # compute the cartesian coordinate and the distance
            distance = struct.dist_r(struct.redCoord[iat], struct.redCoord[jat])

            # xmax act as a cutoff
            if distance > xmax:
                continue

            # select the relevant bond
            for bond in bondsList:
                if struct.atomNames[iat] in bond and struct.atomNames[jat] in bond:
                    selectedBond = bond
                    break

            # add a gaussian
            gaussVal = np.array([g(x, distance, sigma) for x in xval])
            data[selectedBond] += gaussVal

    # print data
    line  = "# structural analysis %s\n" % struct.name
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

if __name__ == "__main__":
    anaStruct()
