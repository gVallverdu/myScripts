#!/usr/bin/env python


# Display settings
import argparse
from pathlib import Path


def add_molecule(cube_file, isovalue=0.02, mode="wireframe", state=1):
    """ Add a molecule in the vmd script """

    if mode == "wireframe":
        mode = 1
        material = "Opaque"
    elif mode == "solid":
        mode = 0
        material = "Transparent"
    else:
        print(f"WARNING: mode '{mode}' is not known, set wireframe")

    new_mol = f"""
mol new {cube_file} type cube first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
# load rep CPK
mol delrep 0 top
mol representation CPK 1.000000 0.500000 20.000000 20.000000
mol color Name
mol selection {{all}}
mol material Opaque
mol addrep top
# load rep Surf pos
mol representation Isosurface {isovalue:.6f} 0 0 {mode} 1 1
mol color ColorID 1
mol selection {{all}}
mol material {material}
mol addrep top
# load rep Surf neg
mol representation Isosurface {-isovalue:.6f} 0 0 {mode} 1 1
mol color ColorID 0
mol selection {{all}}
mol material {material}
mol addrep top
#load rep nodal surface
mol representation Isosurface 0.000000 0 0 1 1 1
mol color ColorID 6
mol selection {{all}}
mol material Opaque
mol addrep top
mol showrep top 3 0
#rename mol
mol rename top {cube_file}
molinfo top set drawn {state}
"""

    return new_mol


def get_options():
    """ Define the argparse parser and returns options """

    def exist(f):
        """ 'Type' for argparse - checks that file exists but does not open """
        if not Path(f).is_file():
            raise argparse.ArgumentTypeError(f"{f} does not exist")
        return f

    # parse command line options
    parser = argparse.ArgumentParser(
        description="VMD Molecular Orbitals",
        epilog="Author: Germain Salvato Vallverdu")
    
    parser.add_argument('-v', '--version', action='version', version="0.1")

    parser.add_argument(
        "-i", "--isovalue", type=float, default=0.02,
        help="Isovalue used to draw the MO."
    )

    group_cube = parser.add_mutually_exclusive_group(required=True)
    group_cube.add_argument(
        "-a", "--all",
        help="Load all cube files in the directory.",
        action="store_const",
        const="all")
    group_cube.add_argument(
        "-c", '--cube',
        help="The cube file containing the OM to be read.",
        metavar="OM.cube",
        type=exist)
    
    group_rep = parser.add_mutually_exclusive_group()
    group_rep.add_argument(
        "-w", "--wireframe",
        help="Draw isosurfaces as wireframes (default).",
        action="store_const",
        const="wireframe")
    group_rep.add_argument(
        "-s", "--solid",
        help="Draw isosurfaces as solid surfaces.",
        action="store_const",
        const="solid")

    parser.add_argument(
        "-o", "--output",
        help="Name of the output VMD script",
        metavar="VMDSCRIPT",
        default="mo.vmd"
    )

    args = parser.parse_args()

    return args


def main():
    """ run the script """

    # get options
    args = get_options()
    
    # vmd file header
    lines = "# VMD script written using save_state $Revision: 1.47 $/n"
    lines += "# VMD version: 1.9.3/n"
    lines += "display projection   Orthographic/n"
    lines += "display depthcue   off/n"
    lines += "display shadows off\n"

    mode = "wireframe" if args.solid is None else "solid"

    if args.all is not None:
        p = Path("./")
        print(f"Read all cub(e) files in: {p.absolute()}")
        cube_files = [f for f in p.glob("*.cub*")]

    else:
        cube_files = [args.cube]

    for i, cube_file in enumerate(cube_files):
        lines += f"# Start with molecule {i}"
        lines += add_molecule(
            cube_file=cube_file, isovalue=args.isovalue, mode=mode,
            state=1 if i == 0 else 0
        )
        lines += f"# done with molecule {i}\n\n"
       
    with open(args.output, "w") as fout:
        fout.write(lines)


if __name__ == '__main__':
    main()
