#!/usr/bin/env python
# -*- coding=utf-8 -*-

"""
BigFiles
--------

print the first bigest files in all subdirectories of a given path

SYNTAX
        bigfiles [OPTIONS] ... [PATH]

DESCRIPTION
        Print the first 20 bigest files in directory [PATH] and in all subdirectories of
        [PATH]. Default [PATH] is the current working directory.

        -h, --help
            print this help

        -n, --nfiles=N
            output the first N bigest files, instead if the first 20. If N is negative,
            print all files.

EXAMPLES
        bigfile -h
        bigfile -n 50
        bigfile ./another/directory/
        bigfile -n -1 
        bigfile -n 20 ./anotherdirectory
"""

__licence__ = "GPL"
__author__ = "Germain Vallverdu <germain.vallverdu@univ-pau.fr>"

import doctest
import os
import sys

KILO_OCTET = 1024
MEGA_OCTET = 1048576
GIGA_OCTET = 1073741824
TERA_OCTET = 1099511627776

def humanSize(taille):
    """ return a byte number in human readable format (1ko 234Mo 2Go

    >>> humanSize( 345)
    '345 o'
    >>> humanSize( 123478)
    '120.6 ko'
    >>> humanSize( 983435743)
    '937.9 Mo'
    >>> humanSize( 12983435743)
    '12.1 Go'
    >>> humanSize( 755812983435743)
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

def BigFiles( dossier = "./", nmax = 20):
    """ print the first bigest files in all subdirectories of a given path 
        
        arguments :
            * dossier : path where you want to list bigest files
            * nmax    : maximum number of files you want BigFiles to list
    """

    fichier = list()
    tailleTotale = 0
    nFiles = 0
    for root, dirs, files in os.walk(dossier):
        for name in files:
            if not os.path.isfile(os.path.join(root,name)): continue
            taille = os.path.getsize(os.path.join(root, name))
            nFiles += 1
            tailleTotale += taille
            fichier.append({"nom":name, "dossier":root, "taille":taille})

    print("\n" + " * * * BigFiles * * *".center(80) + "\n")
    print("--------------------------------------------------------------------------------")
    print(" Analyse du dossier         : " + dossier)
    print(" Nombre de fichiers traités : " + str(nFiles))
    print(" Taille totale des fichiers : " + humanSize(tailleTotale))
    print("--------------------------------------------------------------------------------\n")

    # sort files and print the first nmax
    sorted_fichier = sorted( fichier, key = lambda f:f["taille"], reverse = True)
    for i, f in enumerate(sorted_fichier):
        print(str(i+1).rjust(6) + humanSize(f["taille"]).rjust(10) + f["nom"].rjust(30) + "   " + f["dossier"])
        if i == (nmax-1): break

if __name__ == "__main__":
    doctest.testmod()

    # --------------------------------------------------------
    # default values
    # --------------------------------------------------------
    dossier = "./"
    nmax = 20

    # --------------------------------------------------------
    # get options
    # --------------------------------------------------------
    narg = len(sys.argv)
    if narg > 1:
        args = sys.argv
        i = 1
        while i < narg:
            if args[i] == "-h" or args[i] == "--help":
                print(__doc__)
                exit(0)
            elif args[i] == "-n" or args[i] == "--nfiles":
                try:
                    nmax = int(args[i+1])
                except ValueError:
                    print("\nError : bad arguments " + args[i])
                    print("    try : " + args[0] + " --help\n")
                    exit(1)
                except IndexError:
                    print("\nError : bad arguments " + args[i])
                    print("    try : " + args[0] + " --help\n")
                    exit(1)
                i += 2

            elif args[i] == "-v":
                i += 1
            else:
                if i == narg - 1:
                    dossier = args[i]
                else:
                    print(args)
                    print("\nError : bad arguments")
                    print("    try : " + args[0] + " --help\n")
                    exit(1)
                i += 1

    BigFiles(dossier, nmax)

