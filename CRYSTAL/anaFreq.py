#!/usr/bin/python
# -*- coding=utf-8 -*-

import sys
#from math import sqrt

def readModes(outfile):
    """ read frequencies data on file outfile """

    out = open(outfile, "r")
    line = None
    atom = list()
    while line != "":
        line = out.readline()
        if "NUMBER OF IRREDUCIBLE ATOMS IN THE CONVENTIONAL CELL" in line:
            nat = int(line.split(":")[1])

        if "COORDINATES OF THE EQUIVALENT ATOMS" in line:
            out.readline()
            out.readline()
            out.readline()
            line = out.readline()
            while "NUMBER OF SYMMETRY" not in line:
                atom.append(line.split()[4])
                out.readline()
                line = out.readline()

        if "NORMAL MODES NORMALIZED TO CLASSICAL AMPLITUDES" in line:
            out.readline()
            out.readline()
            out.readline()
            break

    ndep = 3 * nat
    
    nmodes = 0
    modes = list()
    while nmodes != ndep:
        for i in range(6):
            modes.append(list())
        line = out.readline()
        while line.strip() != "":
            if "AT." in line:
                val = [float(v) for i, v in enumerate(line.split()) if i in range(4, 10)]
            else:
                val = [float(v) for v in line.split()[1:]]

            for i in range(6):
                modes[nmodes + i].append(val[i])

            line = out.readline()

        nmodes += 6
        out.readline()
        out.readline()

    out.close()
    print("3 * Natoms      : {0}".format(ndep))
    print("Nombre de modes : {0}".format(nmodes))

    return modes, atom

def readFrequencies(outfile):
    """ read frequencies data on file outfile """

    out = open(outfile, "r")
    line = None
    while line != "":
        line = out.readline()
        if "CONVERSION FACTORS FOR FREQUENCIES:" in line:
            for i in range(7):
                out.readline()
            break

    nmodes = 0
    frequencies = list()
    line = out.readline()
    while line.strip() != "":
        data = [float(val.strip(")")) for i, val in enumerate(line.split()) if i in [2, 3, 4, 9]]
        frequencies.append(data)
        nmodes += 1
        line = out.readline()

    out.close()
    print("Nombre de modes : {0}".format(nmodes))

    # dans frequences : [ua, cm^-1, THz, Intens IR]
    return frequencies

def anaFreq():
    """ build a spectrum """

    outfile = sys.argv[1]
    print("reading file {0}".format(outfile))

    freqData = readFrequencies(outfile)
    modes, atoms = readModes(outfile)

    for n, mode in enumerate(modes):

        tot = 0.
        for val in mode:
            tot += val**2
        #tot = sqrt(tot)

        contrib = dict()
        for i, iat in zip(range(0, 96, 3), range(len(atoms))):
            #at = sqrt(mode[i]**2 + mode[i + 1]**2 + mode[i + 2]**2)
            at = mode[i]**2 + mode[i + 1]**2 + mode[i + 2]**2
            if atoms[iat] not in contrib.keys():
                contrib[atoms[iat]] = at
            else:
                contrib[atoms[iat]] += at

        line = "%12.4f" % freqData[n][1]
        for atom in contrib.keys():
            line += "%10.2f" % (contrib[atom] / tot * 100.0)
        print(line)


        #print("\nFreq %12.4f cm**-1" % freqData[n][1])
    
    for atom in contrib.keys():
        print(atom)

if __name__ == "__main__":
    anaFreq()
