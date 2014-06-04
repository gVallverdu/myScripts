#!/usr/bin/python
# -*- coding=utf-8 -*-

"""
Modules groupes
===============

Contact : Germain Vallverdu <germain.vallverdu@univ-pau.fr>

Description
-----------

Ce module a pour but de decomposer rapidement une representation reductible en
representations irreductibles. Il met a disposition la classe `Groupe` qui contient la
methode `Groupe.reduire(rep)` permettant de reduire une representation et la methode
`Groupe.testGroupe()` qui test les donnees contenu dans l'objet (dimension des RI,
orthognonalite des RI ...).

Tous les groupes ne sont pas disponnibles. Groupes disponnibles :

    * C2v
    * C3v
    * C4v
    * D3h
    * D4h

Exemple d'utilisation
---------------------

>>> from groupes import c4v
>>> rep = [4, 0, 0, 2, 0]
>>> print(c4v)
 C4v  |     1E      2C4      1C2    2sigma_v 2sigma_d
------------------------------------------------------------
  A1  |     1        1        1        1        1    
  A2  |     1        1        1        -1       -1   
  B1  |     1        -1       1        1        -1   
  B2  |     1        -1       1        -1       1    
  E   |     2        0        -2       0        0
>>> c4v.ordre
8
>>> c4v.reduce(rep)
1 A1 + 1 B1 + 1 E

Licence
-------

GPL

"""

import numpy as np

__licence__ = "GPL"
__author__ = "Germain Vallverdu <germain.vallverdu@univ-pau.fr>"
__date__ = "Juin 2012"

class Groupe(object):
    """ groupe ponctuel """

    def __init__(self, name = "", ordre = -1 ):
       """ init group """

       # nom du groupe
       self.name = name

       # ordre du groupe
       self.ordre = ordre

       # Representation irreductible
       self.RI = dict()
       self.ordreRI = list()

       # liste des operations de symetrie
       self.Op = list()

       # classes d'equivalences
       self.classes = dict()
       self.ordreClasses = list()

    def reduce(self, rep):
        """ reduit la representation rep """
        
        if not isinstance(rep, list) and not isinstance(rep, np.ndarray):
            raise TypeError("rep doit etre une liste")

        if len(rep) != len(self.classes.keys()):
            raise ValueError("dimension de rep non coherante")

        nRI = dict()
        for RI in self.RI:
            nRI[RI] = 0
            for i, classe in enumerate(self.ordreClasses):
                nRI[RI] += self.classes[classe] * self.RI[RI][i] * rep[i]
            nRI[RI] /= self.ordre

        ligne = ""
        for RI in self.ordreRI:
            if nRI[RI] == 0:
                continue
            ligne += "{0} {1} + ".format(nRI[RI], RI)
        print(ligne[:-3])

    def testGroupe(self):
        """ Teste les donnees du groupes """

        if self.ordre == -1:
            raise NameError(self.name + " : ordre inconnu")

        if self.name == "":
            raise NameError("nom inconnu")

        if self.classes == {}:
            raise NameError(self.name + " : classes d'equivalence inconnuees")

        if self.ordreClasses == {}:
            raise NameError(self.name + " : ordre des classes d'equivalence inconnuees")

        if self.Op == []:
            raise NameError(self.name + " : Operations inconnues")

        if self.RI == {}:
            raise NameError(self.name + " : RI inconnues")

        if self.ordreRI == []:
            raise NameError(self.name + " : ordre RI inconnues")

        # calcule l'ordre a partir de la definition des classes
        ordre = 0
        for clef in self.classes:
            ordre += self.classes[clef]

        if ordre != self.ordre:
            raise ValueError(self.name + " : Ordre du groupe non coherant avec les classes d'equivalence")

        # verifie l'ordre du groupe et le nombre d'operation
        if ordre != len(self.Op):
            raise ValueError(self.name + " : Ordre du groupe non coherant avec les operations de symetrie")

        # verifie l'ordre des classes
        if len(self.ordreClasses) != len(self.classes.keys()):
            raise ValueError(self.name + " : Nombre de classes incoherent entre l'ordre des classes et les classes")

        # verifie l'ordre des RI
        if len(self.ordreRI) != len(self.RI.keys()):
            raise ValueError(self.name + " : Nombre de RI incoherent entre l'ordre des RI et les RI")

        # verifie les dimensions des RI
        for RI in self.RI:
            if len(self.RI[RI]) != len(self.classes.keys()):
                raise ValueError(self.name + " : RI " + str(RI) + " non coherent avec le nombre de classes")

        # verifie l'orthogonalite des RI
        for RIi in self.ordreRI:
            for RIj in self.ordreRI:
                scal = 0
                for i, classe in enumerate(self.ordreClasses):
                    scal += self.classes[classe] * self.RI[RIi][i] * self.RI[RIj][i]
                if RIi == RIj and not scal / self.ordre == 1:
                    raise ValueError("groupe : " + self.name + " RI " + RIi + " et RI " + RIj + " non normee")
                elif RIi != RIj and scal != 0:
                    raise ValueError("groupe : " + self.name + "RI " + RIi + " et RI " + RIj + " non orthonormee")

        return "ok"

    def __str__(self):
        """ print group """

        ligne = self.name.center(5) + " | "
        for clef in self.ordreClasses:
            ligne += "{0}{1}".format(self.classes[clef], clef).center(9)
        ligne += "\n"

        ligne += "".join(60 * ["-"]) + "\n"

        for RI in self.ordreRI:
            ligne += RI.center(5) + " | "
            for val in self.RI[RI]:
                ligne += str(val).center(9)
            ligne += "\n"

        return ligne

# ----------------------------------------------------------------------------- 
#
# Groupe C2v
#
# ----------------------------------------------------------------------------- 
c2v = Groupe(name = "C2v", ordre = 4)

c2v.classes["E"]        = 1
c2v.classes["C2"]       = 1
c2v.classes["sigma_xz"] = 1
c2v.classes["sigma_yz"] = 1
c2v.ordreClasses = ["E", "C2", "sigma_xz", "sigma_yz"]

c2v.Op = c2v.ordreClasses

c2v.RI["A1"] = [1,  1,  1,  1]
c2v.RI["A2"] = [1,  1, -1, -1]
c2v.RI["B1"] = [1, -1,  1, -1]
c2v.RI["B2"] = [1, -1, -1,  1]

c2v.ordreRI = ["A1", "A2", "B1", "B2"]

# ----------------------------------------------------------------------------- 
#
# Groupe C3v
#
# ----------------------------------------------------------------------------- 
c3v = Groupe(name = "C3v", ordre = 6)

c3v.classes["E"]       = 1
c3v.classes["C3"]      = 2
c3v.classes["sigma_v"] = 3
c3v.ordreClasses = ["E", "C3", "sigma_v"]

c3v.Op = ["E", "C_3^1", "C_3^2", "sigma_1", "sigma_2", "sigma_3"]

c3v.RI["A1"] = [1,  1,  1]
c3v.RI["A2"] = [1,  1, -1]
c3v.RI["E"]  = [2, -1,  0]

c3v.ordreRI = ["A1", "A2", "E"]

# ----------------------------------------------------------------------------- 
#
# Groupe C4v
#
# ----------------------------------------------------------------------------- 
c4v = Groupe(name = "C4v", ordre = 8)

c4v.classes["E"]       = 1
c4v.classes["C4"]      = 2
c4v.classes["C2"]      = 1
c4v.classes["sigma_v"] = 2
c4v.classes["sigma_d"] = 2
c4v.ordreClasses = ["E", "C4", "C2", "sigma_v", "sigma_d"]

c4v.Op = ["E", "C_4^1", "C_4^3", "C2", "sigma_v1", "sigma_v2", "sigma_d1", "sigma_d2"]

c4v.RI["A1"] = [ 1,  1,  1,  1,  1]
c4v.RI["A2"] = [ 1,  1,  1, -1, -1]
c4v.RI["B1"] = [ 1, -1,  1,  1, -1]
c4v.RI["B2"] = [ 1, -1,  1, -1,  1]
c4v.RI["E"]  = [ 2,  0, -2,  0,  0]

c4v.ordreRI = ["A1", "A2", "B1", "B2", "E"]

# ----------------------------------------------------------------------------- 
#
# Groupe D3h
#
# ----------------------------------------------------------------------------- 
d3h = Groupe(name = "D3h", ordre = 12)

d3h.classes["E"]       = 1
d3h.classes["C3"]      = 2
d3h.classes["C2"]      = 3
d3h.classes["sigma_h"] = 1
d3h.classes["S3"]      = 2
d3h.classes["sigma"]   = 3
d3h.ordreClasses = ["E", "C3", "C2", "sigma_h", "S3", "sigma"]

d3h.Op = ["E", "C_3^1", "C_3^2", "C_2a", "C_2a", "C_2a", "sigma_h", "S_3^1", "S_3^5", "sigma_1", "sigma_2", "sigma_3"]

d3h.RI["A1'"]  = [ 1,  1,  1,  1,  1,  1]
d3h.RI["A2'"]  = [ 1,  1, -1,  1,  1, -1]
d3h.RI["E'"]   = [ 2, -1,  0,  2, -1,  0]
d3h.RI["A1''"] = [ 1,  1,  1, -1, -1, -1]
d3h.RI["A2''"] = [ 1,  1, -1, -1, -1,  1]
d3h.RI["E''"]  = [ 2, -1,  0, -2,  1,  0]

d3h.ordreRI = ["A1'", "A2'", "E'", "A1''", "A2''", "E''"]

# ----------------------------------------------------------------------------- 
#
# Groupe D4h
#
# ----------------------------------------------------------------------------- 
d4h = Groupe(name = "D4h", ordre = 16)

d4h.classes["E"]       = 1
d4h.classes["C4"]      = 2
d4h.classes["C2"]      = 1
d4h.classes["C2'"]     = 2
d4h.classes["C2\""]    = 2
d4h.classes["i"]       = 1
d4h.classes["S4"]      = 2
d4h.classes["sigma_h"] = 1
d4h.classes["sigma_v"] = 2
d4h.classes["sigma_d"] = 2
d4h.ordreClasses = ["E", "C4", "C2", "C2'", "C2\"", "i", "S4", "sigma_h", "sigma_v", "sigma_d"]

d4h.Op = ["E", "C_4^1", "C_4^3", "C2", "C2'x", "C2'y", "C2\"", "C2\"", "i", 
          "S_4^1", "S_4^3", "sigma_h", "sigma_v1", "sigma_v2", "sigma_d1", "sigma_d2"]

d4h.RI["A1g"] = [ 1,  1,  1,  1,  1,  1,  1,  1,  1,  1]
d4h.RI["A2g"] = [ 1,  1,  1, -1, -1,  1,  1,  1, -1, -1]
d4h.RI["B1g"] = [ 1, -1,  1,  1, -1,  1, -1,  1,  1, -1]
d4h.RI["B2g"] = [ 1, -1,  1, -1,  1,  1, -1,  1, -1,  1]
d4h.RI["Eg"]  = [ 2,  0, -2,  0,  0,  2,  0, -2,  0,  0]
d4h.RI["A1u"] = [ 1,  1,  1,  1,  1, -1, -1, -1, -1, -1]
d4h.RI["A2u"] = [ 1,  1,  1, -1, -1, -1, -1, -1,  1,  1]
d4h.RI["B1u"] = [ 1, -1,  1,  1, -1, -1,  1, -1, -1,  1]
d4h.RI["B2u"] = [ 1, -1,  1, -1,  1, -1,  1, -1,  1, -1]
d4h.RI["Eu"]  = [ 2,  0, -2,  0,  0, -2,  0,  2,  0,  0]

d4h.ordreRI = ["A1g", "A2g", "B1g", "B2g", "Eg", "A1u", "A2u", "B1u", "B2u", "Eu"]

if __name__ == "__main__":
    """ Test des donnees saisies """

    print(" --- TEST DES DONNEES SAISIEES --- ")
    print("Groupes :")
    print("    {0} : {1}".format("C2v", c2v.testGroupe()))
    print("    {0} : {1}".format("C3v", c3v.testGroupe()))
    print("    {0} : {1}".format("C4v", c4v.testGroupe()))
    print("    {0} : {1}".format("D3h", d3h.testGroupe()))
    print("    {0} : {1}".format("D4h", d4h.testGroupe()))

