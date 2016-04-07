#!/usr/bin/env python
# -*- coding=utf-8 -*-

import os
import sys
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.gridspec import GridSpec

import pymatgen as mg
from pymatgen.io.vasp.outputs import Vasprun, BSVasprun
from pymatgen.electronic_structure.core import Spin, OrbitalType
from pymatgen.electronic_structure.plotter import BSPlotter


class ElStructurePlotter(object):
    """
    Class in order to plot band structure and DOS with an RGB color scale.
    It relies on BSPlotter class of pymatgen.
    """

    def __init__(self, bandpath="./", dospath=None):
        """
        Init method. Read vasprun.xml for the band structure calculation
        and the DOS calculation if a path is provided. The band structure is
        extracted using the fermi level of the dos calculation if available.

        Args:
            bandpath (str): path to vasprun.xml file of the band structure
            dospath (str): path to vasprun.xml file of the dos
        """
        self.xmlbands = os.path.join(bandpath, "vasprun.xml")
        if os.path.exists(self.xmlbands):
            run = BSVasprun(self.xmlbands, parse_projected_eigen=True)
        else:
            raise FileNotFoundError("File {0} not found".format(self.xmlbands))

        kpoints_file = os.path.join(bandpath, "KPOINTS")
        if dospath:
            self.xmldos = os.path.join(dospath, "vasprun.xml")
            if os.path.exists(self.xmldos):
                self.dosrun = Vasprun(self.xmldos)
            else:
                raise FileNotFoundError("File {0} not found".format(self.xmldos))
            self.bands = run.get_band_structure(kpoints_file, line_mode=True,
                                                efermi=self.dosrun.efermi)
        else:
            self.xmldos = None
            self.dosrun = None
            self.bands = run.get_band_structure(kpoints_file, line_mode=True)


def rgbline(ax, k, e, red, green, blue, alpha=1., linewidth=2):
    # creation of segments based on
    # http://nbviewer.ipython.org/urls/raw.github.com/dpsanders/matplotlib-examples/master/colorline.ipynb
    pts = np.array([k, e]).transpose().reshape(-1, 1, 2)
    seg = np.concatenate([pts[:-1], pts[1:]], axis=1)

    nseg = len(k) - 1
    r = [0.5 * (red[i] + red[i + 1]) for i in range(nseg)]
    g = [0.5 * (green[i] + green[i + 1]) for i in range(nseg)]
    b = [0.5 * (blue[i] + blue[i + 1]) for i in range(nseg)]
    a = np.ones(nseg, np.float) * alpha
    lc = LineCollection(seg, colors=list(zip(r, g, b, a)),
                        linewidth=linewidth)
    ax.add_collection(lc)


def make_ticks(ax, bsplot, linewidth=2, linestyle="solid"):
    """
    add ticks and vertical lines to the band structure plot. Add an
    horizontal line at the fermi level.

    Args:
        ax: an Axes object
        bsplot: a BSPlotter object
        linewidth (int): line width
        linestyle (str): line style
    """
    ticks = bsplot.get_ticks()
    uniq_d = list()
    uniq_l = list()
    nticks = len(ticks["label"])
    for i in range(nticks):
        if i == 0:
            uniq_l.append(ticks["label"][i])
            uniq_d.append(ticks["distance"][i])
            ax.axvline(ticks['distance'][i], color='k', linewidth=linewidth,
                       linestyle=linestyle)
        else:
            if ticks["label"][i] == ticks["label"][i - 1]:
                continue
            else:
                uniq_l.append(ticks["label"][i])
                uniq_d.append(ticks["distance"][i])
                ax.axvline(ticks['distance'][i], color='k', linewidth=linewidth,
                           linestyle=linestyle)
    ax.set_xticks(uniq_d)
    ax.set_xticklabels(uniq_l)

    # add fermie level horizontal line
    ax.axhline(y=0, color="k", linewidth=linewidth)


def make_spd_dos_plot(ax, dosrun, spin, reverse=False, **kargs):
    """
    Make a DOS plot

    Args:
        ax: an Axes object
        dosrun: a Vasprun object
        spin: plot Spin.up or Spin.down
        xlim: tuple of plot boundaries
        s, p, d (str): label of projected DOS plot. If None, contribution not plotted
        linewidth (int): line width
        reverse (bool): if True x axes is reversed
    """
    # default options
    s, p, d = None, None, None
    xlim = (0, 5)
    linewidth = 2
    if "s" in kargs:
        s = kargs["s"]
    if "p" in kargs:
        p = kargs["p"]
    if "d" in kargs:
        d = kargs["d"]
    if "xlim" in kargs:
        xlim = kargs["xlim"]
    if "linewidth" in kargs:
        linewidth = kargs["linewidth"]

    # get spd projected DOS
    spd_dos = dosrun.complete_dos.get_spd_dos()

    # title
    if dosrun.is_spin:
        if spin == Spin.up:
            ax.set_title("Spin up", {"fontsize": 18})
        elif spin == Spin.down:
            ax.set_title("Spin down", {"fontsize": 18})
        else:
            raise ValueError("Bad spin value, use Spin.up or Spin.down")

    # tikcs
    if spin == Spin.up:
        ax.set_yticklabels([])
    ax.set_xticklabels([])
    ax.axhline(y=0, color="k", lw=2)
    ax.set_xlabel("Density of States", labelpad=28)
    ax.grid(False)

    # total dos
    ax.fill_between(dosrun.tdos.densities[spin],
                    0,
                    dosrun.tdos.energies - dosrun.efermi,
                    color=(0.7, 0.7, 0.7),
                    facecolor=(0.7, 0.7, 0.7))

    ax.plot(dosrun.tdos.densities[spin],
            dosrun.tdos.energies - dosrun.efermi,
            color=(0.6, 0.6, 0.6),
            label="total DOS")

    # spd contribution
    if s:
        ax.plot(spd_dos[OrbitalType.s].densities[spin],
                dosrun.tdos.energies - dosrun.efermi,
                "r-", label=s, linewidth=linewidth)
    if p:
        ax.plot(spd_dos[OrbitalType.p].densities[spin],
                dosrun.tdos.energies - dosrun.efermi,
                "g-", label=p, linewidth=linewidth)
    if d:
        ax.plot(spd_dos[OrbitalType.d].densities[spin],
                dosrun.tdos.energies - dosrun.efermi,
                "b-", label=d, linewidth=linewidth)

    # reverse x axes
    if reverse:
        ax.set_xlim(xlim[1], xlim[0])
    else:
        ax.set_xlim(xlim[0], xlim[1])


def make_el_dos_plot(ax, dosrun, spin, reverse=False, **kargs):
    """
    Make a DOS plot

    Args:
        ax: an Axes object
        dosrun: a Vasprun object
        spin: plot Spin.up or Spin.down
        xlim: tuple of plot boundaries
        elements: list of pymatgen elements
        linewidth (int): line width
        reverse (bool): if True x axes is reversed
    """
    # default options
    elements = dosrun.initial_structure.composition.elements
    xlim = (0, 5)
    linewidth = 2
    if "elements" in kargs:
        elements = kargs["elements"]
    if "xlim" in kargs:
        xlim = kargs["xlim"]
    if "linewidth" in kargs:
        linewidth = kargs["linewidth"]

    # couleurs  rouge       vert      bleu
    colors = ["#cc0000", "#73d216", "#3465a4"]

    # get spd projected DOS
    eldos = dosrun.complete_dos.get_element_dos()

    # title
    if dosrun.is_spin:
        if spin == Spin.up:
            ax.set_title("Spin up", {"fontsize": 18})
        elif spin == Spin.down:
            ax.set_title("Spin down", {"fontsize": 18})
        else:
            raise ValueError("Bad spin value, use Spin.up or Spin.down")

    # tikcs
    if spin == Spin.up:
        ax.set_yticklabels([])
    ax.set_xticklabels([])
    ax.axhline(y=0, color="k", lw=2)
    ax.set_xlabel("Density of States", labelpad=28)
    ax.grid(False)

    # total dos
    ax.fill_between(dosrun.tdos.densities[spin],
                    0,
                    dosrun.tdos.energies - dosrun.efermi,
                    color=(0.7, 0.7, 0.7),
                    facecolor=(0.7, 0.7, 0.7))

    ax.plot(dosrun.tdos.densities[spin],
            dosrun.tdos.energies - dosrun.efermi,
            color=(0.6, 0.6, 0.6),
            label="total DOS")

    # element dos contribution
    for c, el in zip(colors, elements):
        ax.plot(eldos[el].densities[spin],
                dosrun.tdos.energies - dosrun.efermi,
                color=c,
                linestyle="solid",
                label=el.symbol,
                linewidth=linewidth)

    # reverse x axes
    if reverse:
        ax.set_xlim(xlim[1], xlim[0])
    else:
        ax.set_xlim(xlim[0], xlim[1])


def compute_conrib_spd(bands, spin, name=None):

    # get projected bands
    if name:
        pbands = bands.get_projections_on_elts_and_orbitals({name: ["s", "p", "d"]})
    else:
        raise NameError("Element name is not define")

    # compute s, p, d normalized contributions
    contrib = np.zeros((bands.nb_bands, len(bands.kpoints), 3))
    for b in range(bands.nb_bands):
        for k in range(len(bands.kpoints)):
            sc = pbands[spin][b][k][name]["s"]**2
            pc = pbands[spin][b][k][name]["p"]**2
            dc = pbands[spin][b][k][name]["d"]**2
            tot = sc + pc + dc
            if tot != 0.0:
                contrib[b, k, 0] = sc / tot
                contrib[b, k, 1] = pc / tot
                contrib[b, k, 2] = dc / tot

    return contrib


def compute_contrib_el(bands, spin, elements):

    elbands = bands.get_projection_on_elements()

    # compute normalized contributions
    contrib = np.zeros((bands.nb_bands, len(bands.kpoints), 3))
    for b in range(bands.nb_bands):
        for k in range(len(bands.kpoints)):
            tot = sum([elbands[spin][b][k][el.symbol]**2 for el in elements])
            for i, el in enumerate(elements):
                if tot != 0.0:
                    contrib[b, k, i] = elbands[spin][b][k][el.symbol]**2 / tot

    return contrib


def make_spd_band_plot(ax, bands, yticklabels=True, **kargs):
    """
    Make a DOS plot

    Args:
        ax: an Axes object
        bands: band structure object
        linewidth (int): line width
    """
    # default values of options
    name = None
    linewidth = 2
    if "name"in kargs:
        name = kargs["name"]
    if "linewidth" in kargs:
        linewidth = kargs["linewidth"]
    for key in kargs:
        if key not in ["name", "linewidth"]:
            print("WARNING: option {0} not considered".format(key))

    # band structure plot data
    bsplot = BSPlotter(bands)
    plotdata = bsplot.bs_plot_data(zero_to_efermi=True)

    # spin polarized calculation
    if bands.is_spin_polarized:
        all_spins = [Spin.up, Spin.down]
    else:
        all_spins = [Spin.up]

    for spin in all_spins:
        if spin == Spin.up:
            alpha = 1
            lw = linewidth
        if spin == Spin.down:
            alpha = .7
            lw = linewidth / 2

        # compute s, p, d normalized contributions
        contrib = compute_conrib_spd(bands, spin, name)

        # plot bands
        ikpts = 0
        maxd = -1
        mind = 1e10
        for d, ene in zip(plotdata["distances"], plotdata["energy"]):
            npts = len(d)
            maxd = max(max(d), maxd)
            mind = min(min(d), mind)
            for b in range(bands.nb_bands):
                rgbline(ax, d, ene[str(spin)][b],
                        contrib[b, ikpts:ikpts + npts, 0],
                        contrib[b, ikpts:ikpts + npts, 1],
                        contrib[b, ikpts:ikpts + npts:, 2],
                        alpha, lw)
            ikpts += len(d)

    # add ticks and vlines
    make_ticks(ax, bsplot)
    ax.set_xlabel("k-points")
    ax.set_xlim(mind, maxd)
    ax.grid(False)

    if not yticklabels:
        ax.set_yticklabels([])


def make_el_band_plot(ax, bands, yticklabels=True, **kargs):
    """
    Make a DOS plot

    Args:
        ax: an Axes object
        bands: band structure object
        linewidth (int): line width
    """
    # default values of options
    name = None
    linewidth = 2
    if "name" in kargs:
        name = kargs["name"]
    if "linewidth" in kargs:
        linewidth = kargs["linewidth"]
    if "elements" in kargs:
        elements = kargs["elements"]
    else:
        raise KeyError("argument 'elements' in make_el_band_plot is missing")
    for key in kargs:
        if key not in ["name", "linewidth", "elements"]:
            print("WARNING: option {0} not considered".format(key))

    # band structure plot data
    bsplot = BSPlotter(bands)
    plotdata = bsplot.bs_plot_data(zero_to_efermi=True)

    # spin polarized calculation
    if bands.is_spin_polarized:
        all_spins = [Spin.up, Spin.down]
    else:
        all_spins = [Spin.up]

    for spin in all_spins:
        if spin == Spin.up:
            alpha = 1
            lw = linewidth
        if spin == Spin.down:
            alpha = .7
            lw = linewidth / 2

        # compute s, p, d normalized contributions
        contrib = compute_contrib_el(bands, spin, elements)

        # plot bands
        ikpts = 0
        maxd = -1
        mind = 1e10
        for d, ene in zip(plotdata["distances"], plotdata["energy"]):
            npts = len(d)
            maxd = max(max(d), maxd)
            mind = min(min(d), mind)
            for b in range(bands.nb_bands):
                rgbline(ax, d, ene[str(spin)][b],
                        contrib[b, ikpts:ikpts + npts, 0],
                        contrib[b, ikpts:ikpts + npts, 1],
                        contrib[b, ikpts:ikpts + npts:, 2],
                        alpha, lw)
            ikpts += len(d)

    # add ticks and vlines
    make_ticks(ax, bsplot)
    ax.set_xlabel("k-points")
    ax.set_xlim(mind, maxd)
    ax.grid(False)

    if not yticklabels:
        ax.set_yticklabels([])


def bands_plot(bandpath, dospath, make_bs_plot, make_dos_plot, title=None,
               figsize=(11.69, 8.27), ylim=(-10, 10), bandopt={}, dosopt={},
               legendopt={}):
    """
    Plot the band structure.

    The function assumes a KPOINTS file is present in bandpath folder and the
    band structure calculation was done in line mode.

    Args:
        bandpath (str): path to the vasprun.xml file of the band structure
        dospath (str): path to the vasprun.xml file of the DOS
        make_bs_plot (function): function in order to make the band structure plot
        make_dos_plot (function): function in order to make the DOS plot
        title (str): title of the plot
        figsize (tuple): figure size
        ylim (tuple): y boundaries of the plot
        bandopt (dict): options for the band plot given to make_bs_plot()
        dosopt (dict): options for the dos plot, that are : s, p, d, xlim and linewidth
        legendprop (dict): legend options
    """

    # density of states
    xmldos = os.path.join(dospath, "vasprun.xml")
    print("Reading in file : {0}".format(xmldos))
    dosrun = Vasprun(xmldos)

    # bands
    xmlbands = os.path.join(bandpath, "vasprun.xml")
    print("Reading in file : {0}".format(xmlbands))
    kpoints_file = os.path.join(bandpath, "KPOINTS")
    bandrun = BSVasprun(xmlbands, parse_projected_eigen=True)
    print("Building band structure object")
    bands = bandrun.get_band_structure(kpoints_file,
                                       line_mode=True,
                                       efermi=dosrun.efermi)

    fig = plt.figure(figsize=figsize)
    if not title:
        fig.suptitle("Band structure diagram")
    else:
        fig.suptitle(title)
    if dosrun.is_spin:
        gs = GridSpec(1, 3, width_ratios=[2, 5, 2])
        ax0 = plt.subplot(gs[0])
        ax1 = plt.subplot(gs[1])
        ax2 = plt.subplot(gs[2])
        yticklabels = False
    else:
        gs = GridSpec(1, 2, width_ratios=[2, 1])
        ax1 = plt.subplot(gs[0])
        ax2 = plt.subplot(gs[1])
        yticklabels = True
    gs.update(wspace=0)

    # band structure plot
    print("Making band structure plot")
    make_bs_plot(ax1, bands, yticklabels=yticklabels, **bandopt)

    # Density of states plot
    print("Making DOS plot")
    make_dos_plot(ax2, dosrun, Spin.up, **dosopt)
    if dosrun.is_spin:
        make_dos_plot(ax0, dosrun, Spin.down, **dosopt, reverse=True)

    # plot boundaries and legend
    if dosrun.is_spin:
        ax0.set_ylim(ylim)
        ax0.set_ylabel(r"$E - E_f$   /   eV")
    else:
        ax1.set_ylabel(r"$E - E_f$   /   eV")
    ax1.set_ylim(ylim)
    ax2.set_ylim(ylim)

    ax2.legend(**legendopt)

    print("save fig bsplot.pdf")
    plt.savefig("bsplot.pdf", format="pdf")


if __name__ == "__main__":
    # legend position
    # 'best'         : 0, (only implemented for axes legends)
    # 'upper right'  : 1,
    # 'upper left'   : 2,
    # 'lower left'   : 3,
    # 'lower right'  : 4,
    # 'right'        : 5,
    # 'center left'  : 6,
    # 'center right' : 7,
    # 'lower center' : 8,
    # 'upper center' : 9,
    # 'center'       : 10,

    # set up matplotlib plot
    font = {'family': 'sans', 'size': 20}
    plt.rc('font', **font)
    legendopt = {"fancybox": True, "shadow": True, "prop": {'size': 14},
                 "loc": "center right"}

    # Ni
    root = "../../../Recherche/VASP/TP_VASP/Ni_dev/"
    bandpath = root + "Bandes/"
    dospath = root + "DOS/"
    bandopt = {"name": "Ni"}
    dosopt = {"s": "4s", "p": "4p", "d": "3d", "xlim": (1e-4, 3)}
    bands_plot(bandpath, dospath, make_bs_plot=make_spd_band_plot,
               make_dos_plot=make_spd_dos_plot, legendopt=legendopt,
               bandopt=bandopt, dosopt=dosopt, ylim=(-10, 6))

    # Sr3V2O7
    root = "../../../GitBook/Library/Import/tp-etat-solide/Sr3V2O7/calculs/ucell/"
    bandpath = root + "Bandes_short/"
    dospath = root + "DOS/"
    elements = [mg.Element("Sr"), mg.Element("V"), mg.Element("O")]
    dosopt = {"xlim": (1e-4, 50), "elements": elements}
    bandopt = {"elements": elements}
    bands_plot(bandpath, dospath, make_bs_plot=make_el_band_plot,
               make_dos_plot=make_el_dos_plot, figsize=(11.69, 8.27),
               bandopt=bandopt, dosopt=dosopt, ylim=(-7, 3),
               legendopt=legendopt)
