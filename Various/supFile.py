#!/usr/bin/env python3
# -*- coding=utf-8 -*-

"""
supFile
-------

SYNTAX
    supFile [OPTIONS] 

DESCRIPTION
    Print all files corresponding to the argument list and ask for deletion.
    The script start from the current workind directory and look for all files
    in all subdirectories.
        
    -h, --help
        print this help

EXAMPLES
        supFile CHG WAVECAR
        supFile CHGCAR 
        supFile -h
"""

__licence__ = "GPL"
__author__ = "Germain Vallverdu <germain.vallverdu@univ-pau.fr>"
__date__ = "12 Novembre 2013"

import doctest
import os
import sys
import fnmatch

KILO_OCTET = 1024
MEGA_OCTET = 1048576
GIGA_OCTET = 1073741824
TERA_OCTET = 1099511627776

def humanSize(taille):
    """ return a byte number in human readable format (1ko 234Mo 2Go

    >>> humanSize(345)
    '345 o'
    >>> humanSize(123478)
    '120.6 ko'
    >>> humanSize(983435743)
    '937.9 Mo'
    >>> humanSize(12983435743)
    '12.1 Go'
    >>> humanSize(755812983435743)
    '687.4 To'

    """
    taille = float(taille)
    if taille >= TERA_OCTET:
        size = "%.1f To" % (taille / TERA_OCTET)
    elif taille >= GIGA_OCTET:
        size = "%.1f Go" % (taille / GIGA_OCTET)
    elif taille >= MEGA_OCTET:
        size = "%.1f Mo" % (taille / MEGA_OCTET)
    elif taille >= KILO_OCTET:
        size = "%.1f ko" % (taille / KILO_OCTET)
    else:
        size = "%.0f o" % taille

    return size

def delFile(filenames):
    """ List all file in an path tree and delete it
        
        arguments :
            * filenames (list) : name of the file you want to delete
    """

    # check
    if not isinstance(filenames, list):
        print("filenames must be a list of file name")
        exit(1)

    fichier = list()
    tailleTotale = 0
    nFiles = 0
    for root, dirs, files in os.walk("./"):
        for name in files:
            if not os.path.isfile(os.path.join(root,name)): 
                continue
            for filename in filenames:
                if fnmatch.fnmatch(name, filename):
                    taille = os.path.getsize(os.path.join(root, name))
                    nFiles += 1
                    tailleTotale += taille
                    fichier.append({"nom":name, "dossier":root, "taille":taille})

    print("\n * * * files * * *\n")
    print("--------------------------------------------------------------------------------")
    print(" Nombre de fichiers         : " + str(nFiles))
    print(" Taille totale des fichiers : " + humanSize(tailleTotale))
    print("--------------------------------------------------------------------------------\n")

    # sort files and print 
    sorted_fichier = sorted(fichier, key = lambda f:f["taille"], reverse = True)
    for i, f in enumerate(sorted_fichier):
        print(str(i+1).rjust(6) + humanSize(f["taille"]).rjust(10) + f["nom"].rjust(30) + "   " + f["dossier"])
    print("--------------------------------------------------------------------------------\n")
    print("### WARNING, ALL ABOVE FILES WILL BE DELETE ###")
    answer = ""
    while answer != "y":
        answer = input(" ARE YOU SURE ? (y/n) : ")
        if answer == "n":
            print("nothing was deleted")
            exit(0)
        elif answer == "y":
            continue
        else:
            print(" hit 'y' for yes or 'n' for non")
    for f in sorted_fichier:
        print("rm " + os.path.join(f["dossier"], f["nom"]))
        os.remove(os.path.join(f["dossier"], f["nom"]))
    
if __name__ == "__main__":
    doctest.testmod()

    nargs = len(sys.argv)
    if nargs >= 2:
        if "-h" in sys.argv or "--help" in sys.argv:
            print(__doc__)
            exit(0)
        else:
            filenames = sys.argv[1:]
    else:
        print("Error: You have to give at least one file name")
        print(__doc__)
        exit(1)

    delFile(filenames)

