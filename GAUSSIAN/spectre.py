#!/usr/bin/env python
# -*- coding=utf-8 -*-

"""
Plot an UV-visible spectum from the output of a td-DFT calculation. To each
transition, a gaussian function is added with a given width (in energy) and
an area equal to the oscillator strength.

Syntaxe : 

    spectre.py fichier [sigma, [step ] ]

    fichier : fichier gaussian que l'on va lire (obligatoire)
    sigma   : largeur de la gaussienne (optionnel)
    step    : nombre de points pour construire le spectre (optionnel)

"""
__author__ = "Germain Vallverdu <germain.vallverdu@univ-pau.fr>"
__licence__ = "GPL"
__date__ = "Janvier 2014"

import os
import sys
import scipy as sp
import scipy.constants as cst
import matplotlib.pyplot as plt
from matplotlib.mlab import normpdf

def readTransitions(fichier):
    """ read excitations energies from a gaussian log file 

        :param fichier: Name of the gaussian log file
        :type fichier: string
        :return: List of excitation energies, wavelength and oscillator strength
        :rtype: List
    
    """

    transitions = list()

    td = False
    for line in open(fichier, "r"):
        if "Excitation energies" in line:
            td = True
            continue

        if td:
            if "Excited State" in line:
                ene = float(line.split()[4])
                lam = float(line.split()[6])
                f   = float(line.split("f=")[1])
                transitions.append([ene, lam, f])
            elif "****" in line:
                break

    for i, trans in enumerate(transitions):
        print("%2d   E = %8.4f eV ; L = %8.2f nm ; f = %6.4f" % (i+1, trans[0], trans[1], trans[2]))

    return transitions

def makeSpectre(transitions, sigma, step):
    """ Build a spectrum from transitions energies. For each transitions a gaussian
    function of width sigma is added in order to mimick natural broadening.
    
        :param transitions: list of transitions for readTransitions()
        :type transititions: list
        :param sigma: gaussian width in eV
        :type sigma: float
        :param step: number of absissa value
        :type step: int
        :return: absissa and spectrum value in this order
        :rtype: list, list
    
    """

    # max and min transition energies
    minval = min([val[0] for val in transitions]) - 5.0 * sigma
    maxval = max([val[0] for val in transitions]) + 5.0 * sigma
 
    # points
    npts   = int((maxval - minval) / step) + 1

    # absice
    eneval = sp.linspace(minval, maxval, npts)

    spectre = sp.zeros(npts)
    for trans in transitions:
        spectre += trans[2] * normpdf(eneval, trans[0], sigma)

    return eneval, spectre

def plotSpectre(transitions, eneval, spectre):
    """ plot the UV-visible spectrum using matplotlib. Absissa are converted in nm. """

    # lambda in nm
    lambdaval = [cst.h * cst.c / (val * cst.e) * 1.e9 for val in eneval]

    # plot gaussian spectra
    plt.plot(lambdaval, spectre, "r-", label = "spectre")

    # plot transitions
    plt.vlines([val[1] for val in transitions], \
               0., \
               [val[2] for val in transitions], \
               color = "blue", \
               label = "transitions" )

    plt.xlabel("lambda   /   nm")
    plt.ylabel("Arbitrary unit")
    plt.title("UV-visible spectra")
    plt.grid()
    plt.legend(fancybox = True, shadow = True)
    plt.show()

def spectre(fichier, sigma, step):
    """ call previous routine and make the spectra """

    # read transitions
    transitions = readTransitions(fichier)

    # build spectrum
    eneval, spectre = makeSpectre(transitions, sigma, step)

    # plot spectre
    plotSpectre(transitions, eneval, spectre)

    line  = "# Spectre\n"
    line += "# column 1 : energies in eV\n"
    line += "# column 2 : spectre (arbitrary unit)\n"
    for x, y in zip(eneval, spectre):
        line += "%12.6f %12.6f\n" % (x, y)
    open("spectre.dat", "w").write(line)
    
if __name__ == "__main__":

    fichier = sys.argv[1]
    if not os.path.exists(fichier):
        print("File %s does not exist" % fichier)

    try:
        sigma   = float(sys.argv[2])
    except IndexError:
        sigma = 0.05

    try:
        step = float(sys.argv[3])
    except IndexError:
        step = 0.001

    print("--------------------------------------")
    print("File  : %s" % fichier)
    print("sigma : %8.4f" % sigma)
    print("step  : %8.4f" % step)
    print("--------------------------------------")

    spectre(fichier, sigma, step)

