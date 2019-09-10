#!/usr/bin/env python
# -*- coding=utf-8 -*-

"""
 _   _   _____   _____   __   _   _____   _____
| | | | |  _  \ /  _  \ |  \ | | |  _  \ /  ___/
| |_| | | |_| | | | | | |   \| | | | | | | |___
|  _  | |  _  { | | | | | |\   | | | | | \___  \\
| | | | | |_| | | |_| | | | \  | | |_| |  ___| |
|_| |_| |_____/ \_____/ |_|  \_| |_____/ /_____/


Hydrogen bond analysis from a GROMACS trajectory. Use the following command
in order to get a formated and centered GROMACS trajectory from the traj_comp.xtc
file.

    gmx trjconv -pbc mol -center -ur compact -s run.tpr -f traj_comp.xtc -o trajout.gro

The list of hydrogen bonds you want to investigate have to be defined in a file,
for example `hb.yml`. An sample file is given as example.

"""

# Dumping schematic of time series after each h-bond, key follows:
#    |          .       -       o      x      *      @    |
#       0-5%   5-20%  20-40%  40-60% 60-80% 80-95% 95-100% occupancy

import re
import yaml
import argparse
import time
import io
import mdtraj as md

import numpy as np

TITLE = """
 _   _   _____   _____   __   _   _____   _____
| | | | |  _  \ /  _  \ |  \ | | |  _  \ /  ___/
| |_| | | |_| | | | | | |   \| | | | | | | |___
|  _  | |  _  { | | | | | |\   | | | | | \___  \\
| | | | | |_| | | |_| | | | \  | | |_| |  ___| |
|_| |_| |_____/ \_____/ |_|  \_| |_____/ /_____/
"""


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

def compute_hbonds(trajectory, data, donors, acceptors, solvent_donors, solvent_acceptors):
    """ compute hbonds """

    # parameters for HBonds calculations
    rcut = data["rcut"]
    acut = np.radians(data["angle"])

    # number of frame, for array dimention
    nframe = trajectory.xyz.shape[0]

    # dict to sotre time series
    hbond_data = dict()

    # case of solvent
    if "solvent" in data:
        solvent_data = dict()
        for iat in donors:
            solvent_data[(iat,)] = np.full((nframe, 3), np.nan)
            solvent_data[(iat,)][:, 0].fill(0.)
        for acceptor in acceptors:
            iat = acceptor[0]
            for ih in acceptor[1:]:
                solvent_data[(iat, ih)] = np.full((nframe, 3), np.nan)
                solvent_data[(iat, ih)][:, 0].fill(0.)

    print("\nStart computing hydrogen bonds")
    print("------------------------------")

    t1 = time.time()
    t0 = t1
    for iframe, coords in enumerate(trajectory.xyz):

        if iframe % (nframe // 10) == 0 and iframe != 0:
            t2 = time.time()
            dt = t2 - t1
            print(" **** %3.0f%% in %.0f s - frame %4d (ETA %.0f s)" % (
                iframe / nframe * 100, dt, iframe, dt * 10 - (t2 - t0)))
            t1 = t2

        for donor in donors:
            for acceptor in acceptors:

                # manage atom who are both acceptor and donor
                if donor == acceptor[0]:
                    continue

                rij = coords[donor] - coords[acceptor[0]]
                dij = np.sqrt((rij**2).sum())

                if dij < rcut:
                    for ih in acceptor[1:]:
                        rik = coords[ih] - coords[acceptor[0]]
                        dik = np.sqrt((rik**2).sum())
                        angle = np.arccos(np.dot(rij, rik) / dij / dik)

                        if angle < acut:
                            # hbond exist
                            if (donor, (acceptor[0], ih)) in hbond_data:
                                hbond_data[(donor, (acceptor[0], ih))][iframe, 0] += 1
                                hbond_data[(donor, (acceptor[0], ih))][iframe, 1] = dij
                                hbond_data[(donor, (acceptor[0], ih))][iframe, 2] = angle
                            else:
                                hbond_data[(donor, (acceptor[0], ih))] = np.full(
                                    (nframe, 3), np.nan)
                                hbond_data[(donor, (acceptor[0], ih))][:, 0].fill(0.)
                                hbond_data[(donor, (acceptor[0], ih))][iframe, 0] = 1
                                hbond_data[(donor, (acceptor[0], ih))][iframe, 1] = dij
                                hbond_data[(donor, (acceptor[0], ih))][iframe, 2] = angle

        if "solvent" in data:
            # solvent acceptor
            for s_acceptor in solvent_acceptors:
                for donor in donors:
                    if donor == s_acceptor[0]:
                        continue

                    rij = coords[donor] - coords[s_acceptor[0]]
                    dij = np.sqrt((rij**2).sum())

                    if dij < rcut:
                        # compute angle for all H
                        for ih in s_acceptor[1:]:
                            rik = coords[ih] - coords[s_acceptor[0]]
                            dik = np.sqrt((rik**2).sum())
                            angle = np.arccos(np.dot(rij, rik) / dij / dik)

                            if angle < acut:
                                # hbond exist
                                solvent_data[(donor,)][iframe, 0] += 1
                                solvent_data[(donor,)][iframe, 1] = dij
                                solvent_data[(donor,)][iframe, 2] = angle

            # solvent donor
            for s_donor in solvent_donors:
                for acceptor in acceptors:
                    if s_donor == acceptor[0]:
                        continue

                    rij = coords[s_donor] - coords[acceptor[0]]
                    dij = np.sqrt((rij**2).sum())

                    if dij < rcut:
                        # compute angle
                        for ih in acceptor[1:]:
                            rik = coords[ih] - coords[acceptor[0]]
                            dik = np.sqrt((rik**2).sum())
                            angle = np.arccos(np.dot(rij, rik) / dij / dik)

                            if angle < acut:
                                # hbond exist
                                solvent_data[(acceptor[0], ih)][iframe, 0] += 1
                                solvent_data[(acceptor[0], ih)][iframe, 1] = dij
                                solvent_data[(acceptor[0], ih)][iframe, 2] = angle

    dt = time.time() - t0
    print("Done in %.0f s" % dt)

    if "solvent" in data:
        return hbond_data, solvent_data
    else:
        return hbond_data, None

def get_symbol(val):
    """  Returns the symbol corresponding to the given value.
       |          .       -       o      x      *      @    |
          0-5%   5-20%  20-40%  40-60% 60-80% 80-95% 95-100%

        Args:
            val (float): a number in [0, 1]

        Returns:
            symbol (str)
    """
    if not 0 <= val <= 1:
        raise ValueError("value is not in [0, 1]")

    symbols = [" ", ".", "-", "o", "x", "*", "@"]
    step = 1 / 5

    i = int(val / step) + 1

    if val < 0.05:
        i = 0
    elif val >=0.95:
        i = 6

    return symbols[i]

def print_res(topology, hbond_data):
    """ Output a table with the results """

    print("\nRESULTS")
    print("-------")

    print("\nDumping schematic of time series after each h-bond, key follows:")
    print("   |          .       -       o      x      *      @    |")
    print("      0-5%   5-20%  20-40%  40-60% 60-80% 80-95% 95-100% occupancy")

    nframe = len(hbond_data[list(hbond_data.keys())[0]][:, 0])
    window = nframe // 10 # windows in order to output time series with symbols

    print("\n| ------------------------------------------------------------------------------------------ |")
    print("|      DONOR     |    ACCEPTORH       ACCEPTOR    |                                          |")
    print("| atom# :res@atom| atom# :res@atom atom# :res@atom| %occupied   distance         angle       |")
    print("| ------------------------------------------------------------------------------------------ |")
    line_ts = "# time series for hbonds\n"
    line_ts += "# column  1: step\n"
    icol = 1
    all_time_series = [np.arange(nframe)]
    for donor, (acceptor, acceptorh) in hbond_data:
        a_donor = topology.atom(donor)
        a_acceptor = topology.atom(acceptor)
        h_acceptor = topology.atom(acceptorh)

        icol += 1
        line_ts += "# column %2d: " % icol

        col = "%5d :%s@%s" % (a_donor.serial, str(a_donor.residue), a_donor.name)
        line = "|" + col.ljust(16)
        line_ts += col + " --> "
        col = "%5d :%s@%s" % (h_acceptor.serial, str(h_acceptor.residue), h_acceptor.name)
        line += "|" + col.ljust(16)
        line_ts += col
        col = "%5d :%s@%s" % (a_acceptor.serial, str(a_acceptor.residue), a_acceptor.name)
        line += col.ljust(16) + "|"
        line_ts += col + "\n"

        occupency = hbond_data[donor, (acceptor, acceptorh)][:, 0].mean() * 100
        distance = hbond_data[donor, (acceptor, acceptorh)][:, 1]
        angle = np.degrees(hbond_data[donor, (acceptor, acceptorh)][:, 2])
        line += " %6.2f %8.3f (%5.2f) %8.3f (%5.2f)" % (occupency,
                                                        np.nanmean(distance), np.nanstd(distance),
                                                        np.nanmean(angle), np.nanstd(angle))

        line += " |"

        time_serie = hbond_data[donor, (acceptor, acceptorh)][:, 0]
        all_time_series.append(time_serie)

        # output average time series by window
        time_serie_window = [time_serie[i * window: (i+1) * window].mean() for i in range(9)]
        time_serie_window.append(time_serie[8 * window:].mean())
        line += "".join([get_symbol(tsw) for tsw in time_serie_window])

        line += "|"
        print(line)
    print("| ------------------------------------------------------------------------------------------ |")

    # save time series to file
    str_time_series = io.BytesIO()
    all_time_series = np.array(all_time_series).transpose()
    np.savetxt(str_time_series, all_time_series, delimiter=" ", fmt="%5d")

    line_ts += str_time_series.getvalue().decode()
    with open("hbond_time_series.dat", "w") as f:
        f.write(line_ts)

def print_solvent(topology, s_data, data):
    """ Output a table with the results about hydrogen bonds with solvent """

    print("\n| ------------------------------------------------------------------------------------------ |")
    print("|  Hydrogen bonds with solvent molecules")
    if "donors" in data["solvent"]:
        print("|  Solvent donors:")
        for donor in data["solvent"]["donors"]:
            print("|    - %s" % donor)
    if "acceptors" in data["solvent"]:
        print("|  Solvent acceptors:")
        for acceptor in data["solvent"]["acceptors"]:
            print("|    - %s %s %s" % tuple(acceptor))


    print("| ------------------------------------------------------------------------------------------ |")
    print("|      DONOR     |    ACCEPTORH       ACCEPTOR    |                                          |")
    print("| atom# :res@atom| atom# :res@atom atom# :res@atom|  ave nbr    distance         angle       |")
    print("| ------------------------------------------------------------------------------------------ |")
    for iat in s_data:
        if len(iat) == 1:
            atom = topology.atom(iat[0])

            col = "%5d :%s@%s" % (iat[0], str(atom.residue), atom.name)
            line = "|" + col.ljust(16)
            line += "|" + 32 * " " + "|"

            occupency = s_data[iat][:, 0].mean()
            distance = s_data[iat][:, 1]
            angle = np.degrees(s_data[iat][:, 2])
            line += " %6.2f %8.3f (%5.2f) %8.3f (%5.2f)" % (occupency,
                                                            np.nanmean(distance), np.nanstd(distance),
                                                            np.nanmean(angle), np.nanstd(angle))
        else:
            a_atom = topology.atom(iat[0])
            h_atom = topology.atom(iat[1])
            line = "|" + 16 * " "
            col = "%5d :%s@%s" % (iat[0], str(a_atom.residue), a_atom.name)
            line += "|" + col.ljust(16)
            col = "%5d :%s@%s" % (iat[1], str(h_atom.residue), h_atom.name)
            line += col.ljust(16) + "|"

            occupency = s_data[iat][:, 0].mean()
            distance = s_data[iat][:, 1]
            angle = np.degrees(s_data[iat][:, 2])
            line += " %6.2f %8.3f (%5.2f) %8.3f (%5.2f)" % (occupency,
                                                            np.nanmean(distance), np.nanstd(distance),
                                                            np.nanmean(angle), np.nanstd(angle))


        line += " |"
        print(line)
    print("| ------------------------------------------------------------------------------------------ |")

def init(data, topology):
    """ Set up parameters and print info """

    # print parameters
    print("\nCut-off")
    print("-------")
    print("Rcut:  %5.2f nm" % data["rcut"])
    print("Angle: %5.2f (deg)" % data["angle"])

    # print list of donors and acceptors
    if "acceptors" in data:
        print("\nList of acceptors:")
        print("------------------")
        for acceptor in data["acceptors"]:
            ha = topology.atom(acceptor[0] - 1)
            line = "%4s(%4d) @ %4s(%4d) :" % (ha.name, ha.serial, ha.residue.name,
                                        ha.residue.index + 1)
            for ih in acceptor[1:]:
                h = topology.atom(ih - 1)
                line += " %s(%d)" % (h.name, h.serial)
            print(line)
        # !!!! switch to index [0 -> N]
        acceptors = [[iat - 1 for iat in acceptor] for acceptor in data["acceptors"]]

    else:
        acceptors = []
        print("WARNING: you did not define any acceptor.")

    if "donors" in data:
        print("\nList of donors:")
        print("----------------")
        for donor in data["donors"]:
            hd = topology.atom(donor - 1)
            line = "%4s(%4d) @ %4s(%4d)" % (hd.name, hd.serial, hd.residue.name,
                    hd.residue.index + 1)
            print(line)
         # !!!! switch to index [0 -> N]
        donors = [donor - 1 for donor in data["donors"]]
    else:
        donors = []
        print("WARNING: you did not define any donor.")

    if donors == [] and acceptors == []:
        raise KeyError("You did not define neither donor nor acceptor.")

    # set up in case of solvent
    if "solvent" in data:
        soldata = data["solvent"]

        if "name" in soldata:
            solres = soldata["name"]
        else:
            raise KeyError("You did not give the name of solvent residues.")

        if "donors" in soldata:
            sol_donor_names = soldata["donors"]
        else:
            sol_donor_names = []
            print("Warning: You did not define any solvent donors")

        if "acceptors" in soldata:
            sol_acceptor_names = soldata["acceptors"]
        else:
            sol_acceptor_names = []
            print("Warning: You did not define any solvent acceptors")

        if sol_donor_names == [] and sol_acceptor_names == []:
            raise KeyError("You did not define neither donors nor acceptors for the solvent.")

        # list of solvent donors
        solvent_donors = list()
        for atom in topology.atoms:
            if atom.residue.name == solres and atom.name in sol_donor_names:
                solvent_donors.append(atom.index)
        solvent_donors = np.array(solvent_donors)

        # list of solvent acceptors
        solvent_acceptors = list()
        for res in topology.residues:
            if res.name == solres:
                for sol_acceptor in sol_acceptor_names:
                    ih = list()
                    for atom in res.atoms:
                        if atom.name == sol_acceptor[0]:
                            iacceptor = [atom.index]
                        if atom.name in sol_acceptor[1:]:
                            ih.append(atom.index)
                solvent_acceptors.append(iacceptor + ih)
        solvent_acceptors = np.array(solvent_acceptors)

    else:
        solvent_donors = []
        solvent_acceptors = []

    return donors, acceptors, solvent_donors, solvent_acceptors

def get_options():
    """ get options from the command line """

    parser = argparse.ArgumentParser(
        prog="hbonds",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=__doc__
    )

    # trying to mimic gromacs gmx option names
    parser.add_argument("-i", "--input",
                        help="input file name describing hbond pairs (default hb.yml)",
                        metavar="INPUT.yml",
                        default="hb.yml",
                        type=str)

    parser.add_argument("-f", "--trajfile",
                        help="GROMACS trajectory file.",
                        metavar="TRAJFILE",
                        default="trajout.gro",
                        type=str)

    parser.add_argument("-s", "--top",
                        help="A file which contains topology: a gro or pdf file for example.",
                        metavar="TOP",
                        default="confout.gro",
                        type=str)

    parser.add_argument("-r", "--rcut",
                        help="cut off distance for hydrogen bonds in nm (default 0.32 nm)",
                        metavar="R (nm)",
                        default=0.32,
                        type=float)

    parser.add_argument("-a", "--angle",
                        help="cut off angle for hydrogen bonds in degrees (default 40°)",
                        metavar="ALPHA",
                        default=40,
                        type=float)

    return vars(parser.parse_args())


if __name__ == "__main__":

    print(TITLE)

    # read command line arguments
    args = get_options()
    with open(args["input"], "r") as f:
        data = yaml.load(f)
    data.update(args)

    # read trajectory
    trajectory = read_xtc(**args)

    # set up data
    donors, acceptors, solvent_donors, solvent_acceptors = init(data, trajectory.top)

    # compute hbonds
    hbond_data, solvent_data = compute_hbonds(trajectory, data, donors, acceptors,
                                              solvent_donors, solvent_acceptors)

    # print results
    print_res(trajectory.top, hbond_data)
    if "solvent" in data:
        print_solvent(trajectory.top, solvent_data, data)
