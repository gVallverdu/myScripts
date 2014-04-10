#!/usr/bin/env python
# -*- coding=utf-8 -*-

import sys
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.gridspec import GridSpec

import pymatgen as mg
from pymatgen.symmetry.bandstructure import HighSymmKpath
from pymatgen.io.vaspio.vasp_output import Vasprun
from pymatgen.electronic_structure.core import Spin

def rgbline(ax, k, e, red, green, blue, alpha=1.):
    #creation of segments based on
    #http://nbviewer.ipython.org/urls/raw.github.com/dpsanders/matplotlib-examples/master/colorline.ipynb
    pts = np.array([k, e]).T.reshape(-1, 1, 2)
    seg = np.concatenate([pts[:-1], pts[1:]], axis=1)

    nseg = len(k)-1
    r = [0.5*(red[i]+red[i+1]) for i in range(nseg)]
    g = [0.5*(green[i]+green[i+1]) for i in range(nseg)]
    b = [0.5*(blue[i]+blue[i+1]) for i in range(nseg)]
    a = np.ones(nseg, np.float)*alpha
    lc = LineCollection(seg, colors=zip(r,g,b,a), linewidth = 2)
    ax.add_collection(lc)

if __name__ == "__main__":
    # read data
    # ---------

    # kpoints labels
    path = HighSymmKpath(mg.read_structure("./POSCAR")).kpath["path"]
    labels = list()
    for p in path:
        if len(labels) != 0:
            labels[-1] = r"$%s,%s$" % (labels[-1].strip("$"), p[0])
        else:
            labels.append(r"$%s$" % p[0])
        for kp in p[1:]:
            labels.append(r"$%s$" % kp)

    # density of state
    dosrun = Vasprun("../../DOS/Cu/vasprun.xml")
    spd_dos = dosrun.complete_dos.get_spd_dos()

    # bands
    run = Vasprun("vasprun.xml", parse_projected_eigen = True)
    bands = run.get_band_structure("KPOINTS", 
                                   line_mode = True, 
                                   efermi = dosrun.efermi)

    # set up matplotlib plot
    # ----------------------

    # general options for plot
    font = {'family': 'serif', 'size': 24}
    plt.rc('font', **font)

    # set up 2 graph with aspec ration 2/1
    # plot 1: bands diagram
    # plot 2: Density of State
    gs = GridSpec(1, 2, width_ratios=[2,1])
    fig = plt.figure(figsize=(11.69, 8.27))
    fig.suptitle("Bands diagram of copper")
    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1]) # , sharey=ax1)

    # set ylim for the plot
    # ---------------------
    emin = -10.
    emax = 10.
    ax1.set_ylim(emin, emax)
    ax2.set_ylim(emin, emax)

    # Band Diagram
    # ------------
    name = "Cu"
    pbands = bands.get_projections_on_elts_and_orbitals({name: ["s", "p", "d"]})

    # compute s, p, d normalized contributions
    contrib = np.zeros((bands.nb_bands, len(bands.kpoints), 3))
    for b in range(bands.nb_bands):
        for k in range(len(bands.kpoints)):
            sc = pbands[Spin.up][b][k][name]["s"]**2
            pc = pbands[Spin.up][b][k][name]["p"]**2
            dc = pbands[Spin.up][b][k][name]["d"]**2
            tot = sc + pc + dc
            if tot != 0.0:
                contrib[b, k, 0] = sc / tot
                contrib[b, k, 1] = pc / tot
                contrib[b, k, 2] = dc / tot
            
    # plot bands using rgb mapping
    for b in range(bands.nb_bands):
        rgbline(ax1, 
                range(len(bands.kpoints)),
                [e - bands.efermi for e in bands.bands[Spin.up][b]],
                contrib[b,:,0], 
                contrib[b,:,1], 
                contrib[b,:,2])

    # style
    ax1.set_xlabel("k-points")
    ax1.set_ylabel(r"$E - E_f$   /   eV")
    ax1.grid()

    # fermi level at 0
    ax1.hlines(y=0, xmin=0, xmax=len(bands.kpoints), color="k", lw=2)

    # labels
    nlabs = len(labels)
    step = len(bands.kpoints) / (nlabs - 1)
    for i, lab in enumerate(labels):
        ax1.vlines(i * step, emin, emax, "k")
    ax1.set_xticks([i * step for i in range(nlabs)])
    ax1.set_xticklabels(labels)

    ax1.set_xlim(0, len(bands.kpoints))

    # Density of state
    # ----------------

    ax2.set_yticklabels([])
    ax2.grid()
    ax2.set_xlim(1e-6, 32)
    ax2.hlines(y=0, xmin=0, xmax=32, color="k", lw=2)
    ax2.set_xlabel("Density of State")

    # spd contribution
    ax2.plot(spd_dos["S"].densities[Spin.up],
        dosrun.tdos.energies - dosrun.efermi, 
        "r-", label = "3s", lw = 2)
    ax2.plot(spd_dos["P"].densities[Spin.up], 
        dosrun.tdos.energies - dosrun.efermi, 
        "g-", label = "3p", lw = 2)
    ax2.plot(spd_dos["D"].densities[Spin.up],
        dosrun.tdos.energies - dosrun.efermi, 
        "b-", label = "3d", lw = 2)

    # total dos
    ax2.fill_between(dosrun.tdos.densities[Spin.up],  
        0, 
        dosrun.tdos.energies - dosrun.efermi, 
        color = (0.7, 0.7, 0.7),
        facecolor = (0.7, 0.7, 0.7))

    ax2.plot(dosrun.tdos.densities[Spin.up], 
        dosrun.tdos.energies - dosrun.efermi,
        color = (0.6, 0.6, 0.6),
        label = "total DOS")

    # plot format style
    # -----------------
    ax2.legend(fancybox=True, shadow=True, prop={'size':18})
    plt.subplots_adjust(wspace = 0)

    #plt.show()
    plt.savefig(sys.argv[0].strip(".py") + ".pdf", format="pdf")

