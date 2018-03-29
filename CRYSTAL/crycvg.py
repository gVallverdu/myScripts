#!/usr/bin/env python3
# -*- coding=utf-8 -*-

"""
crycvg

Read a CRYSTAL output file and plot convergence data.
"""

__author__ = "Germain Vallverdu"
__email__ = "germain.vallverdu@univ-pau.fr"
__licence__ = "GPL"

import sys
import matplotlib.pyplot as plt
from crystalio import CrystalOutfile


def crycvg(filename):
    """
    Plot relevant quantities in order to follow a geometry convergence with CRYSTAL
    """

    cryout = CrystalOutfile(filename)
    nstep = len(cryout.convergence_data)

    print("title             : ", cryout.title)
    print("Optimization steps: ", nstep)

    plt.rcParams.update({"axes.grid": True, "font.size": 16})

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2, sharex=True,
                                                 figsize=(12, 8))
    fig.suptitle("CRYSTAL: Convergence panel")

    ax1.plot([e - min(cryout.energies) for e in cryout.energies], marker="o",
             label="$E_{min}$ = %e ua" % min(cryout.energies), linestyle="--")
    ax1.set_ylabel("$E - E_{min}$ (ua)")
    ax1.set_xlabel("Optimization steps")
    ax1.legend()

    ax2.plot([data["max_grad"] for data in cryout.convergence_data],
             marker="o", linestyle="dashed", label="max gradient")
    ax2.plot([data["rms_grad"] for data in cryout.convergence_data],
             label="rms gradient", marker="o", linestyle="dashed")
    ax2.plot([data["max_grad_thr"] for data in cryout.convergence_data], color="C0")
    ax2.plot([data["rms_grad_thr"] for data in cryout.convergence_data], color="C1")
    ax2.set_ylabel("MAX/RMS GRADIENT")
    ax2.set_xlabel("Optimization steps")
    ax2.legend()

    # displacement start at step 1
    ax3.plot(range(1, nstep), [data["max_displac"] for data in cryout.convergence_data[1:]],
             label="Max Displacement", marker="o", linestyle="dashed")
    ax3.plot(range(1, nstep), [data["rms_displac"] for data in cryout.convergence_data[1:]],
             label="RMS Displacement", marker="o", linestyle="dashed")
    ax3.plot(range(1, nstep), [data["max_displac_thr"] for data in cryout.convergence_data[1:]],
             color="C0")
    ax3.plot(range(1, nstep), [data["rms_displac_thr"] for data in cryout.convergence_data[1:]],
             color="C1")
    ax3.set_ylabel("MAX/RMS Displacement")
    ax3.set_xlabel("Optimization steps")
    ax3.legend()

    # Gradient norm is not printed if STEP is rejected !
    plt_data = list()
    for step, data in enumerate(cryout.convergence_data):
        if "norm_grad" in data:
            plt_data.append([step, data["norm_grad"], data["norm_grad_thr"]])
        else:
            plt.axvline(step, color="C3")
    ax4.plot([data[0] for data in plt_data], [data[1] for data in plt_data],
             label="norm grad", marker="o", linestyle="dashed")
    ax4.plot([data[0] for data in plt_data], [data[2] for data in plt_data], color="C0")
    ax4.set_ylabel("NORME GRADIENT")
    ax4.set_xlabel("Optimization steps")
    ax4.text(1, 1, "Rejected steps", color="C3", horizontalalignment='right',
             verticalalignment='bottom', transform=ax4.transAxes, fontsize=12)

    # fig.tight_layout()
    # ax1.set_xticklabels(ax1.get_xticks())
    # [tk.set_visible(True) for tk in ax1.get_xticklabels()]

    plt.show()


if __name__ == "__main__":
    crycvg(sys.argv[1])
