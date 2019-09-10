#!/usr/bin/env python
# -*- coding=utf-8 -*-

"""
                    ___
                   |__ \\
   __ _  _ __  ___    ) | __ _  _ __ ___
  / _  ||  __|/ _ \  / / / _  ||  _   _ \\
 | (_| || |  | (_) |/ /_| (_| || | | | | |
  \__, ||_|   \___/|____|\__, ||_| |_| |_|
   __/ |                    | |
  |___/                     |_|

This program export a Gromacs configuration from a centered trajectory to an
input file for a QM calculation. Use the following command in order to get a
centered GROMACS trajectory from the traj_comp.xtc file.

    gmx trjconv -pbc mol -center -ur compact -s run.tpr -f traj_comp.xtc -o trajout.xtc


"""

import os
import re
import argparse
import time

import numpy as np

import pymatgen as mg
from pymatgen.io.gaussian import GaussianInput

import mdtraj as md

TITLE = """
                    ___
                   |__ \\
   __ _  _ __  ___    ) | __ _  _ __ ___
  / _  ||  __|/ _ \  / / / _  ||  _   _ \\
 | (_| || |  | (_) |/ /_| (_| || | | | | |
  \__, ||_|   \___/|____|\__, ||_| |_| |_|
   __/ |                    | |
  |___/                     |_|

"""


def extract_confs(trajectory, **kwargs):
    """
    Extract a snapshot and make a QM input file. In order to set up the QM input
    file you have to write the custum_work function. That function must take
    as arguments a numpy array of the coordinates of the current snapshot and a
    mdtraj.Topology object. Moreover, the whole input arg are passed as a kwargs.
    """

    # number of frame, for array dimention
    nframe, natom = trajectory.xyz.shape[:2]

    # get needed variables:
    interval = kwargs["interval"]
    path = kwargs["path"]
    begin = kwargs["begin"]
    end = kwargs["end"]
    qm_software = kwargs["qm_software"]

    print("Start configuration extraction")
    print("------------------------------")

    t1 = time.time()
    t0 = t1
    jobs = ""
    for iframe, coords in enumerate(trajectory.xyz):

        # exclude frame before 'begin'
        if iframe < begin:
            continue

        # exclude frame after 'end'
        if end != -1 and iframe > end:
            continue

        # time evaluation
        if iframe % (nframe // 10) == 0 and iframe != 0:
            t2 = time.time()
            dt = t2 - t1
            print(" **** %3.0f%% in %.0f s - frame %4d (ETA %.0f s)" % (
                iframe / nframe * 100, dt, iframe, dt * 10 - (t2 - t0)))
            t1 = t2

        # extract one frame every interval
        if (iframe - begin) % interval == 0:

            qm_coords, qm_names = custom_work(coords, mdtop=trajectory.top, **kwargs)

            # make the input file
            if qm_software == "gaussian":
                jobs += write_gau_input(qm_names, qm_coords, iframe, **kwargs)
            elif qm_software == "orca":
                print("Will be supported soon !")
                raise ValueError()
            else:
                raise ValueError("Wrong QM software name.")

            jobs += "\n"

    dt = time.time() - t0
    print("Done in %.0f s" % dt)

    # script to launch jobs
    with open(os.path.join(path, "lancer.sh"), "w") as f:
        f.write("#!/bin/bash\n\n")
        f.write("cd %s\n\n" % os.path.join(os.getcwd(), path))
        f.write(jobs)


def write_gau_input(species, coords, iframe, jobname="run_gau", **kwargs):
    """ write a gaussian input file """

    path = kwargs["path"]
    functional = kwargs["functional"]
    basis_set = kwargs["basis_set"]
    nprocs = kwargs["nprocs"]

    # get first letter of atom name => element
    species = [specie[0] for specie in species]

    # convert nm to angstrom
    coords = np.array(coords) * 10

    # build Molecule and Gaussian object
    mol = mg.Molecule(species, coords, validate_proximity=True)
    ginp = GaussianInput(
        mol,
        charge=0,
        title="Written by gro2qm - snapshot %d" % iframe,
        functional=functional, basis_set=basis_set,
        route_parameters={"pop": "(MK, MBS)"},
        link0_parameters={"%nprocs": nprocs},
        dieze_tag="#"
    )

    basename = "gau_%05d" % iframe
    ginp.write_file(os.path.join(path, basename + ".com"), cart_coords=True)
    mol.to(fmt="xyz", filename=os.path.join(path, basename + ".xyz"))

    return "sbatch -J snap_%05d %s %s.com" % (iframe, jobname, basename)


def read_xtc(trajfile, top, **kwargs):
    """ read the trajectrory using mdtraj module from the xtc file """

    print("Read trajectory file: ", trajfile)
    print("Read topology from:   ", top)

    t1 = time.time()
    trajectory = md.load(trajfile, top=top)
    t2 = time.time() - t1

    nframes, natoms = trajectory.xyz.shape[:2]
    print("# frames: ", nframes)
    print("# atoms:  ", natoms)
    print("Read trajecory in %.0f s\n" % t2)

    return trajectory


def get_options():
    """ get options from the command line """

    parser = argparse.ArgumentParser(
        prog="gro2gau",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=__doc__
    )

    # trying to mimic gromacs gmx option names
    parser.add_argument("-f", "--trajfile",
                        help="GROMACS trajectory file.",
                        metavar="TRAJFILE",
                        default="traj_comp.xtc",
                        type=str)

    parser.add_argument("-s", "--top",
                        help="A file which contains topology: a gro or pdb file for example.",
                        metavar="TOP",
                        default="confout.gro",
                        type=str)

    parser.add_argument("-p", "--path",
                        help="output path for QM input files",
                        metavar="PATH",
                        default="./gau/",
                        type=str)

    parser.add_argument("-i", "--interval",
                        help="Frame is extracted if mod(iframe, interval) == 0",
                        metavar="I",
                        default=10,
                        type=int)

    parser.add_argument("-b", "--begin",
                        help="Index of the first snapshot.",
                        metavar="I",
                        default=0,
                        type=int)

    parser.add_argument("-e", "--end",
                        help="Index of the last snapshot.",
                        metavar="I",
                        default=-1,
                        type=int)

    # QM calculation options
    parser.add_argument("-q", "--qm_software",
                        help="Name of the QM software for which input files will be done (gaussian or orca).",
                        metavar="QM_SOFTWARE",
                        default="gaussian",
                        choices=["gaussian", "orca"],
                        type=str)

    parser.add_argument("--functional",
                        help="Hamiltonian for the QM calculation.",
                        metavar="FUNC",
                        default="cam-B3LYP",
                        type=str)

    parser.add_argument("--basis_set",
                        help="Basis set for the QM calculation.",
                        metavar="BASIS",
                        default="6-31+G**",
                        type=str)

    parser.add_argument("--nprocs",
                        help="Number of cpu required for the QM calculation.",
                        metavar="PROCS",
                        default=1,
                        type=int)

    return vars(parser.parse_args())


def custom_work(coords, mdtop, **kwargs):
    """ """

    # data
    rcut = 0.7
    natmol1 = 1 - 1
    natmol2 = 26 - 1

    # conformation C
    #  2 = N2
    #  9 = O3
    # 12 = O4
    # 16 = O1
    # 20 = O2
    # 24 = N1
    # hetatoms = [i - 1 for i in [2, 9, 12, 16, 20, 24]]
    # nwat_by_atoms = [3, 1, 3, 3, 3, 3]

    # conformation B
    #  2 = N2
    #  9 = O3
    # 12 = O1
    # 16 = O2
    # 25 = O4
    # 20 = N1
    hetatoms = [i - 1 for i in [2, 9, 12, 16, 25, 20]]
    nwat_by_atoms = [3, 1, 3, 3, 3, 3]

    # coordinates of the molecule
    qm_coords = coords[natmol1: natmol2 + 1]
    qm_names = [mdtop.atom(iat).name for iat in range(natmol1, natmol2 + 1)]
    qm_coords = qm_coords.tolist()

    # look for the first sovation shell of the molecule
    dmax = 0.
    nwat = 0
    all_water_selected = list()
    for hetatom, nwat_by_atom in zip(hetatoms, nwat_by_atoms):
        water_selected = list()
        for atom in mdtop.atoms:
            # exclude the solute molecule
            if natmol1 <= atom.index <= natmol2:
                continue

            # look for O atom of water molecule
            if atom.name == "O" and atom.residue.name == "HOH" and \
               atom.index not in all_water_selected:

                iat = atom.index
                rij = coords[iat] - coords[hetatom]
                dij = np.sqrt((rij**2).sum())

                if dij <= rcut:
                    water_selected.append((iat, dij))

        # extract the nwat_by_atom cloesest molecules
        water_selected.sort(key=lambda item: item[1])
        d_val = ""
        nLH = 0
        for iwat in range(nwat_by_atom):
            iat, d = water_selected[iwat]

            # add to the whole water selected
            all_water_selected.append(iat)

            nwat += 1
            dmax = max(dmax, d)
            if d < 0.32:
                d_val += "%8.3f*" % d
                nLH += 1
            else:
                d_val += "%8.3f " % d

            # add the water molecule to qm_coords
            qm_coords.append(coords[iat].tolist())
            qm_coords.append(coords[iat + 1].tolist())
            qm_coords.append(coords[iat + 2].tolist())
            qm_names += ["OW", "HW", "HW"]

        print("{:4d} {:4d} {:4d} {}".format(hetatom, nwat_by_atom, nLH, d_val))


    print("{:10.4f} {:5d}".format(dmax, nwat))
    return qm_coords, qm_names

if __name__ == "__main__":

    # read command line arguments
    args = get_options()
    if not os.path.exists(args["path"]):
        os.mkdir(args["path"])
        print("Make folder: ", args["path"])
    if not os.path.isdir(args["path"]):
        raise OSError("%s is not a directory." % args["path"])

    print(TITLE)

    # read trajectory
    trajectory = read_xtc(**args)

    # make gaussian inputs
    extract_confs(trajectory, **args)

# add = False
# if add:
#     # add H
#     # N2, H14, H13, C6
#     #  1,   0,   2,  3 => isomere C
#     N = 1
#     H1 = 0
#     H2 = 2
#     C = 3
#     b1 = coords[H1, :] - coords[N, :]
#     b1 /= np.sqrt((b1**2).sum())
#     b2 = coords[H2, :] - coords[N, :]
#     b2 /= np.sqrt((b2**2).sum())
#     b3 = coords[C, :] - coords[N, :]
#     b3 /= np.sqrt((b3**2).sum())
#     vec = -(b1 + b2 + b3)
#     vec /= np.sqrt((vec**2).sum())
#     newH = coords[N, :] + 1 * vec
#
#     species.append("H")
#     coords = np.concatenate([coords, newH.reshape(1, 3)], axis=0)
