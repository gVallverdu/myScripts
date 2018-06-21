#!/usr/bin/env python

"""
opt_slab\n\n

A quick script to write VASP inputs to explore the lattice parameters of a
slab system.
"""

import argparse
import shutil
from pathlib import Path
import numpy as np

# use pymatgen objects
import pymatgen as mg
from pymatgen.io.vasp.inputs import Poscar


def copy_files(wd, args, struct):
    """ Copy files into the give directory """

    shutil.copy(args.incar, wd)
    shutil.copy(args.kpoints, wd)
    shutil.copy(args.potcar, wd)
    shutil.copy(args.jobfile, wd / ("j" + wd.name))

    Poscar(struct).write_file(str(wd / "POSCAR"))


def scale_structure(struct, sf):
    """
    Scale the lattice of a structure according to the choosen dimension and
    using the give scale factor.
    """

    sf_array = np.array(sf)

    # multiply by line
    matrix = struct.lattice.matrix
    matrix = (matrix.T * sf_array).T

    return mg.Structure(mg.Lattice(matrix), struct.species, struct.frac_coords)


def gen_calcs(args):
    """ set up all calculations """

    # set up scale factors and dimension
    if args.from_file:
        with open(args.from_file, "r") as f:
            scale_factors = list()
            dim = list()
            for line in f:
                values = [float(val) for val in line.split()]
                scale_factors.append(values)
                dim.append(len(values))

        if not all(d == dim[0] for d in dim):
            raise ValueError("You have to provide the same number of scaling "
                             "factors on each line.")
        else:
            dimension = args.dimension
            if dim[0] == 1:
                if dimension > 1 and dimension != -3:
                    raise ValueError("scale factors values and dimension are not consistent.")
            elif dim[0] != dimension:
                raise ValueError("scale factors values and dimension are not consistent.")

            for sf in scale_factors:
                if dimension == 1:
                    sfx = sf[0]
                    if args.axis == "a":
                        sf += [1.0, 1.0]
                    elif args.axis == "b":
                        sf.append(1.0)
                        sf.insert(0, 1.0)
                    elif args.axis == "c":
                        sf.insert(0, 1.0)
                        sf.insert(0, 1.0)
                elif dimension == 2:
                    sf += [1.]
                elif dimension == -3:
                    sfx = sf[0]
                    sf += [sfx, sfx]
            scale_factors = [tuple(sf) for sf in scale_factors]

    else:
        dimension = args.dimension
        scale_factors = list()
        for sfx in args.scale_factors:
            if dimension == 1:
                if args.axis == "a":
                    scale_factors.append((sfx, 1.0, 1.0))
                elif args.axis == "b":
                    scale_factors.append((1.0, sfx, 1.0))
                elif args.axis == "c":
                    scale_factors.append((1.0, 1.0, sfx))
                continue
            elif dimension == -3:
                scale_factors.append((sfx, sfx, sfx))
                continue
            for sfy in args.scale_factors:
                if dimension == 2:
                    scale_factors.append((sfx, sfy, 1.0))
                    continue
                for sfz in args.scale_factors:
                    scale_factors.append((sfx, sfy, sfz))

    # set up working directory name
    if dimension == 1:
        wd = Path(args.to) / "scan_axes"
    elif dimension == 2:
        wd = Path(args.to) / "scan_surf"
    elif dimension == 3:
        wd = Path(args.to) / "scan_anisotrop"
    elif dimension == -3:
        wd = Path(args.to) / "scan_vol"

    wd.mkdir(parents=True)

    print("Files will be generated in:", str(wd))
    print(" * Initial structure:", args.poscar)
    print(" * INCAR file from:", args.incar)
    print(" * KPOINTS file from:", args.kpoints)
    print(" * POTCAR file from:", args.potcar)
    print(" * job file:", args.jobfile)

    # initial structure
    struct0 = mg.Structure.from_file(args.poscar)

    # bash script to submit all Files
    lancer = "#!/bin/bash\n\n"
    scale_values = "# scaling factor values and folder names\n"

    # make inputs
    print("\nGrid in %d dimension(s)" % dimension)
    print("Number of calculations: ", len(scale_factors))
    print(scale_factors)

    for sfx, sfy, sfz in scale_factors:

        if dimension == 1 or dimension == -3:
            val = sfx if args.axis == "a" else sfy if args.axis == "b" else sfz
            wd_sf = wd / ("SF_%.2f" % val)
            print("--> Run %5.2f in %s" % (val, wd_sf))
            scale_values += "%10f %s\n" % (val, wd_sf.name)
        elif dimension == 2:
            wd_sf = wd / ("SF_%.2f_%.2f" % (sfx, sfy))
            print("--> Run (%.2f, %.2f) in %s" % (sfx, sfy, wd_sf))
            scale_values += "%10f %10f %s\n" % (sfx, sfy, wd_sf.name)
        elif dimension == 3:
            wd_sf = wd / ("SF_%.2f_%.2f_%.2f" % (sfx, sfy, sfz))
            print("--> Run (%.2f, %.2f, %.2f) in %s" % (sfx, sfy, sfy, wd_sf))
            scale_values += "%10f %10f %10f %s\n" % (sfx, sfy, sfz, wd_sf.name)

        wd_sf.mkdir()

        structi = scale_structure(struct0, (sfx, sfy, sfz))
        copy_files(wd_sf, args, structi)

        lancer += "cd %s\n" % str(wd_sf.resolve())
        lancer += "sbatch j%s\n\n" % wd_sf.name

    with open(wd / "lancer.sh", "w") as f:
        f.write(lancer)

    with open(wd / "scale_values.dat", "w") as f:
        f.write(scale_values)


def get_options():
    """ get options from command lines """

    parser = argparse.ArgumentParser(prog="opt_slab", description=__doc__)

    # scaling factor to be used
    parser.add_argument(
        "-s", "--scale_factors",
        type=float, nargs="+",
        default=[0.99, 1.0, 1.01],
        help="Scaling factors of lattice parameters",
        metavar="SCALE_FACTOR"
    )

    # dimension
    parser.add_argument(
        "-d", "--dimension", type=int, default=1,
        help="Dimension of the grid 1, 2 or 3.\n If dimension is 2, be sure the"
             " slab is perpendicular to the c axis. Dimension 3 is for an "
             "anisotropic grid. Dimension -3 is for an isotropic scan.",
        choices=(1, 2, 3, -3), metavar="DIMENSION"
    )

    # dimension
    parser.add_argument(
        "-a", "--axis", type=str, default="a",
        help="If dimension is 1, define the axis along which the scan is done.",
        choices=("a", "b", "c"), metavar="AXIS"
    )

    # scaling factors from a file
    parser.add_argument(
        "-f", "--from_file", type=str, default=None, metavar="FILE",
        help="The list of scaling factors is provided in a file. "
             "Each line corresponds to one calculation and provides 1, 2 or 3 "
             "values depending on the dimension.",
    )

    # workind directory for all calculations
    parser.add_argument("-t", "--to",
                        help="parent directory where inputs will be written",
                        metavar="WORKDIR",
                        default="./",
                        type=str)

    # vasp input files
    parser.add_argument("-i", "--incar", metavar="INCAR",
                        default="./INCAR", type=str,
                        help="INCAR file to be used for the calculations.")

    parser.add_argument("-p", "--poscar", metavar="POSCAR",
                        default="./POSCAR", type=str,
                        help="POSCAR file to be used for the calculations.")

    parser.add_argument("-k", "--kpoints", metavar="KPOINTS",
                        default="./KPOINTS", type=str,
                        help="KPOINTS file to be used for the calculations.")

    parser.add_argument("-c", "--potcar", metavar="POTCAR",
                        default="./POTCAR", type=str,
                        help="POTCAR file to be used for the calculations.")

    parser.add_argument("-j", "--jobfile", metavar="JOBFILE",
                        default="./job", type=str,
                        help="job file to be used in order to submit the calculations.")

    return parser.parse_args()


if __name__ == "__main__":
    args = get_options()
    gen_calcs(args)
