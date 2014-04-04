#!/usr/bin/env python
# -*- coding=utf-8 -*-

""" lecture de l'output de CRYSTAL et construction de la surface et du POSCAR Pour VASP """

__author__ = "Germain Vallverdu"
__email__ = "germain.vallverdu@univ-pau.fr"
__licence__ = "GPL"

import crystal
import sys

# --------------------
# elements
# --------------------
elementsName = ["H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al",
"Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni",
"Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc",
"Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La"]

# --------------------
# hauteur de vide en A
# --------------------
zvide = 18.0

# ------------------------
# lecture output CRYSTAL09
# ------------------------
lignes = open(sys.argv[1], "r").readlines()

sectionCoord = False
tmp = list()
for i, ligne in enumerate(lignes):

    if "DEFINITION OF THE NEW LATTICE VECTORS" in ligne:
        while "LATTICE PARAMETERS  (ANGSTROM" not in ligne:
            i += 1
            ligne = lignes[i]
        valeurs = [float(val) for val in lignes[i+2].split() ]
        slab = crystal.Crystal(a = valeurs[0], \
                               b = valeurs[1], \
                               c = valeurs[2], \
                               alpha = valeurs[3], \
                               beta  = valeurs[4], \
                               gamma = valeurs[5], \
                               name = "slab CRYSTAL")
        print(slab)

    if "COORDINATES OF THE ATOMS BELONGING TO THE SLAB" in ligne:
        sectionCoord = True
        continue

    if sectionCoord:
        # on saute la premiere ligne
        if "LAB" in ligne:
            continue

        # si ligne blanche on sort
        if ligne.strip() == "" :
            break

        z = int(ligne.split()[1])
        valeurs = [ float(val) for val in ligne[56:].split() ]
        tmp.append( [z, valeurs[0], valeurs[1], valeurs[2] ] )

# --------------------
# cherche zmin et zmax
# --------------------
zmin = 1.e15
zmax = -1.e15
for x in tmp:
    if x[3] < zmin:
        zmin = x[3]
    if x[3] > zmax:
        zmax = x[3]
print("zmin  = {:12.7f}".format(zmin))
print("zmax  = {:12.7f}".format(zmax))
print("dz    = {:12.7f}".format(zmax - zmin))

# ---------------------------------
# choix de la nouvelle valeurs de c
# ---------------------------------
newc = zmax - zmin + zvide
print("new c = {:12.7f}\n".format(newc))

# -----------------------------------
# nouveau slab avec c perpendiculaire
# -----------------------------------

# nombres d'atomes
slab.Natoms = len(tmp)
print("nombre atomes lus = " + str(slab.Natoms) )

newSlab = crystal.Crystal(a = slab.a, \
                          b = slab.b, \
                          c = newc, \
                          alpha = 90., \
                          beta  = 90.,  \
                          gamma = slab.gamma, \
                          name = "new slab for VASP")
newSlab.Natoms = len(tmp)
print(newSlab)

# --------------------------------------------
# decalage des coordonnees entre 0 et 1
# conversion de z en coordonnees reduites
# decalage de la surface au centre de la boite
# --------------------------------------------
for iat in range(newSlab.Natoms):
    tmp[iat][1] += 0.5
    tmp[iat][2] += 0.5
    tmp[iat][3] = tmp[iat][3] / newSlab.c + 0.5

# ----------------------
# composition du slab
# ----------------------
atoms = dict()
for iat in range(newSlab.Natoms):
    if elementsName[tmp[iat][0] - 1] not in atoms.keys():
        atoms[elementsName[tmp[iat][0] - 1]] = 1
    else:
        atoms[elementsName[tmp[iat][0] - 1]] += 1
    
for element in atoms.keys():
    print("Nbre de {:2s} = {:d}".format(element, atoms[element]))
    for i in range(atoms[element]):
        newSlab.atomNames.append(element)

# ----------------------------------
# on replace les atomes dans l'ordre
# ----------------------------------
redCoord = list()
nat = 0
for element in atoms.keys():
    z = elementsName.index(element) + 1
    for iat in range(newSlab.Natoms):
        if tmp[iat][0] == z:
            redCoord.append(tmp[iat][1:])
            nat += 1

if nat != len(tmp):
    print("Error atom number")
    exit(1)

newSlab.redCoord = redCoord

# ------------------
# Ecriture du POSCAR 
# ------------------
newSlab.toPOSCAR()

