#!/usr/bin/env python
# -*-coding:utf-8 -*-

"""
Apply scofield cross sections to a DOS
Journal of Electron Spectroscopy and Related Phenomena, 8 (1976) 129-137)
"""

from __future__ import division, print_function

import numpy as np
import matplotlib.pyplot as plt
from pymatgen.electronic_structure.core import Spin
from pymatgen.io.vasp.outputs import Vasprun

__author__ = "Germain Salvato-Vallverdu"
__email__ = "germain.vallverdu@univ-pau.fr"
__licence__ = "GPL"
__date__ = "June 2014"

def modulate(xmlfile="vasprun.xml", tofile=False):
    """
    Multiply DOS by cross sections and plot it using matplotlib

    Args:
        xmlfile(string): name of the vasprun.xml file
        tofile(bool): data are printed into file BV.dat if true
    """

    run = Vasprun(xmlfile)
    dos = run.complete_dos

    BV = dict()
    BV["total"] = np.zeros(dos.energies.shape)
    for el in dos.structure.composition.elements:
        # cross section
        sigma = CrossSec.from_string(el.symbol)
        print(sigma.comment)

        # spd DOS of the element
        el_dos = dos.get_element_spd_dos(el)

        # sum up spin contribution
        spd_dos = dict()
        if run.is_spin:
            for orb in el_dos.keys():
                spd_dos[orb] = el_dos[orb].densities[Spin.up] \
                             + el_dos[orb].densities[Spin.down]
        else:
            for orb in el_dos.keys():
                spd_dos[orb] = el_dos[orb].densities[Spin.up]

        # compute the BV
        BV[el] = sigma.s * spd_dos["S"] \
               + sigma.p * spd_dos["P"] \
               + sigma.d * spd_dos["D"]

        BV["total"] += BV[el]

    if tofile:
        lines = "# Valence Band of compound %s\n" % \
            dos.structure.composition.reduced_formula
        lines += "# column 1: E - E_fermi (eV)\n"
        lines += "# column 2: total valence band\n"
        i = 2
        for el in dos.structure.composition.elements:
            i += 1
            lines += "# column %d: Contribution of %s \n" % (i, el.symbol)
        for i in xrange(len(dos.energies)):
            lines += "%12.7f " % (dos.energies[i] - dos.efermi)
            lines += "%12.7f " % BV["total"][i]
            for el in dos.structure.composition.elements:
                lines += "%12.7f " % BV[el][i]
            lines += "\n"

        open("BV.dat", "w").write(lines)

    else:
        # plot
        font = {'family': 'serif', 'size': 20}
        plt.rc('font', **font)
        plt.figure(figsize=(11.69, 8.27)) # A4

        plt.fill_between(dos.energies - dos.efermi, BV["total"], 0,
            color=(.8, .8, .8))
        plt.plot(dos.energies - dos.efermi, BV["total"], color=(.5, .5, .5),
            label="total")
        color = ["r-", "g-", "b-", "y-", "m-", "c-"]
        for el, c in zip(dos.structure.composition.elements, color):
            plt.plot(dos.energies - dos.efermi, BV[el], c, label=el.symbol)

        # set up figure
        plt.title("Valence band of " + \
            run.initial_structure.composition.reduced_formula)
        plt.xlabel(r"$E - E_{f}$   /   eV")
        plt.ylabel("Density of states times cross sections")
        plt.grid()
        ymin, ymax = plt.ylim()
        plt.vlines(0, ymin, ymax, color="k", lw=1, linestyle="--")
        plt.legend(prop={"size":18})

        plt.show()

class Kalpha(object):
    """
    Enum type for energies of K alpha photons
    """

    Mg = 1254
    Al = 1487
    all_energies = {"Mg": Mg, "Al": Al}

class CrossSecSpecie(object):
    """
    Store cross sections values of a specie
    """

    def __init__(self, ener, name, s=0, p1=0, p2=0, d1=0, d2=0, comment = None):
        """
        Store cross sections values of element of symbol name

        Args:
            ener: Photon energy
            name: Element symbol
            s: cross section for s electrons
            p1: cross section for p1/2 electrons
            p2: cross section for p3/2 electrons
            d1: cross section for d3/2 electrons
            d2: cross section for d5/2 electrons
        """

        if ener not in Kalpha.all_energies.values():
            print("Error, energy must be in ", Kalpha.all_energies)
            exit(1)

        self._ener = ener
        self._name = name
        self._s = s
        self._p1 = p1
        self._p2 = p2
        self._d1 = d1
        self._d2 = d2
        if not comment:
            comment = name + " cross section for energy " + str(ener)
        self._comment = comment

    @property
    def energy(self):
        """
        Photon energy
        """
        return self._ener

    @property
    def name(self):
        """
        specie name
        """
        return self._name

    @property
    def s(self):
        """
        p1/2 cross section
        """
        return self._s

    @property
    def p_half(self):
        """
        p1/2 cross section
        """
        return self._p1

    @property
    def p_three_halves(self):
        """
        p3/2 cross section
        """
        return self._p2

    @property
    def d_three_halves(self):
        """
        d3/2 cross section
        """
        return self._d1

    @property
    def d_five_halves(self):
        """
        d5/2 cross section
        """
        return self._d2

    @property
    def p(self):
        """
        Cross section for p1/2 and p3/2 electrons weighted by the multiplicity.
        """
        return (2*self._p1 + 4*self._p2)/6

    @property
    def d(self):
        """
        Cross section for d3/2 and d5/2 electrons weighted by the multiplicity.
        """
        return (4*self._d1 + 6*self._d2)/10

    @property
    def comment(self):
        """
        Comment on value
        """
        return self._comment

    def __repr__(self):
       return self.name

    def __str__(self):
       line = "Specie name: " + self._name + "\n"
       line += "Photon energy: " + str(self._ener) + "\n"
       line += self._comment + "\n"
       line += "s = " + str(self.s) + "\n"
       line += "p = " + str(self.p) + "\n"
       line += "d = " + str(self.d) + "\n"
       return line

class CrossSec(object):
    """
    Cross section object from publication of J.H. Scofield (Journal of Electron
    Spectroscopy and Related Phenomena, 8 (1976) 129-137).
    """

    # second row
    Li_p_Mg = CrossSecSpecie(Kalpha.Mg, "Li_p_Mg", 0.0593,
        comment="Li+, 1s2 configuration, Mg Kalpha")
    Li_p_Al = CrossSecSpecie(Kalpha.Al, "Li_p_Al", 0.0568,
        comment="Li+, 1s2 configuration, Al Kalpha")

    N_Mg = CrossSecSpecie(Kalpha.Mg, "N_Mg", 0.0841, 0.0025, 0.0049,
        comment="N: 1s2 2s2 2p3 configuration, Mg Kalpha")
    N_Al = CrossSecSpecie(Kalpha.Al, "N_Al", 0.0867, 0.0022, 0.0043,
        comment="N: 1s2 2s2 2p3 configuration, Al Kalpha")

    O_Mg = CrossSecSpecie(Kalpha.Mg, "O_Mg", 0.1345, 0.0073, 0.0145,
        comment="O: 1s2 2s2 2p4 configuration, Mg Kalpha")
    O_Al = CrossSecSpecie(Kalpha.Al, "O_Al", 0.1405, 0.0065, 0.0128,
        comment="O: 1s2 2s2 2p4 configuration, Al Kalpha")

    # third row
    Na_p_Mg = CrossSecSpecie(Kalpha.Mg, "Na_p_Mg", 0.390, 0.0714, 0.1406,
        comment="Na+: 1s2 2s2 2p6 configuration, Mg Kalpha")
    Na_p_Al = CrossSecSpecie(Kalpha.Al, "Na_p_Al", 0.422, 0.0654, 0.1287,
        comment="Na+: 1s2 2s2 2p6 configuration, Al Kalpha")

    P_Mg = CrossSecSpecie(Kalpha.Mg, "P_Mg", 0.0998, 0.0129, 0.0253,
        comment="P: 1s2 2s2 2p6 3s2 3p3 configuration, Mg Kalpha")
    P_Al = CrossSecSpecie(Kalpha.Al, "P_Al", 0.1116, 0.0124, 0.0244,
        comment="P: 1s2 2s2 2p6 3s2 3p3 configuration, Al Kalpha")

    S_Mg = CrossSecSpecie(Kalpha.Mg, "S_Mg", 0.1302, 0.0269, 0.0527,
        comment="S: 1s2 2s2 2p6 3s2 3p4 configuration, Mg Kalpha")
    S_Al = CrossSecSpecie(Kalpha.Al, "S_Al", 0.1465, 0.0262, 0.0512,
        comment="S: 1s2 2s2 2p6 3s2 3p3 configuration, Al Kalpha")

    # fourth row
    Mn_Mg = CrossSecSpecie(Kalpha.Mg, "Mn_Mg", 0.0398, 0, 0, 0.0484, 0.0711,
        comment="Mn: 4s2 3d5, Mg Kalpha, 0 for p cros section")
    Mn_Al = CrossSecSpecie(Kalpha.Al, "Mn_Al", 0.0464, 0, 0, 0.0424, 0.0622,
        comment="Mn: 4s2 3d5, Al Kalpha, 0 for p cros section")
    Fe_Mg = CrossSecSpecie(Kalpha.Mg, "Fe_Mg", 0.0425, 0, 0, 0.0788, 0.1156,
        comment="Cu: 4s2 3d6, Mg Kalpha, 0 for p cros section")
    Fe_Al = CrossSecSpecie(Kalpha.Al, "Fe_Al", 0.0497, 0, 0, 0.1017, 0.0497,
        comment="Cu: 4s2 3d6, Al Kalpha, 0 for p cros section")
    Cu_Mg = CrossSecSpecie(Kalpha.Mg, "Cu_Mg", 0.0188, 0, 0, 0.268, 0.390,
        comment="Cu: 4s1 3d10, Mg Kalpha, 0 for p cros section")
    Cu_Al = CrossSecSpecie(Kalpha.Al, "Cu_Al", 0.0221, 0, 0, 0.240, 0.349,
        comment="Cu: 4s1 3d10, Al Kalpha, 0 for p cros section")

    all_data = (Li_p_Mg, Li_p_Al, N_Mg, N_Al, O_Mg, O_Al,
                Na_p_Mg, Na_p_Al, P_Mg, P_Al, S_Mg, S_Al,
                Mn_Mg, Mn_Al, Fe_Mg, Fe_Al, Cu_Mg, Cu_Al)

    @staticmethod
    def from_string(specie, ener=Kalpha.Al):
        """
        Return cross section of a specie from a string representation.
        """
        for cs in CrossSec.all_data:
            if specie in cs.name and ener == cs.energy:
                return cs
        raise ValueError("%s not found in the data base" % specie)

if __name__ == "__main__":
    xmlfile = "vasprun.xml"
    modulate(xmlfile, tofile=True)
