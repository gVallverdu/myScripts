# coding: utf-8

"""
This script provide a way to parse GROMACS hydrogen bond data from an XPM file,
an index file, and a GROMACS .gro file. It outputs a summary of hydrogen bonds
including donor and acceptor atom names, indices, and the number of hydrogen
bonds observed in each frame.

WARNING: Be careful about the ordering of hydrogen bonds in the XPM file.
The aim of the xpm file is to output a picture of the hydrogen bonds. Thus
the y-axis is increasing from bottom to top. In consequence, the first line
of the xpm file corresponds to the last bond entry in the index file.

https://gromacs.org-gmx-users.maillist.sys.kth.narkive.com/PX932N3W/gmx-users-ordering-of-hydrogen-bonds-in-hbn-and-hbm-output-in-g-hbond
"""

import sys

from pathlib import Path
import numpy as np


def parse_xpm_hbond_matrix(xpm_file: str):
    """Parses a GROMACS XPM file containing hydrogen bond data and returns a binary matrix."""

    meta_data = False

    p = Path(xpm_file)
    if not p.is_file():
        raise FileNotFoundError(f"File {xpm_file} does not exist.")

    with p.open("r") as f:
        for line in f:
            if line.startswith("/*"):
                continue
            if line.startswith('"'):
                try:
                    vals = [int(val) for val in line.strip('",\n').split()]
                    meta_data = True
                except TypeError as error:
                    print(f"Error parsing line: {line}. Error: {error}")
                break

        if not meta_data:
            return matrix

        # skip lines of color definitions
        [f.readline() for _ in range(vals[2])]

        lines = list()
        for line in f:
            if line.startswith("/*"):
                continue
            else:
                line = line.strip('",\n')

                if len(line) == vals[0]:
                    lines.append(line)
                else:
                    raise ValueError(
                        f"Line length mismatch: expected {vals[0]}, got {len(line)}"
                    )

    matrix = np.array([list(line) for line in lines])
    matrix = (matrix == "o").astype(np.uint8)
    return matrix


def read_gro(gro_file: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Reads a GROMACS .gro file and returns the coordinates, atom names
    and box dimensions.
    In the context of this script the aim is only to get atom names.

    Returns:
        tuple: A tuple containing:
            - atom_names (np.ndarray): Array of atom names.
            - coords (np.ndarray): Array of coordinates (x, y, z).
            - box (np.ndarray): Box dimensions as a 3-element array.
    """
    gro_path = Path(gro_file)
    if not gro_path.is_file():
        raise FileNotFoundError(f"File {gro_file} does not exist.")

    with gro_path.open("r") as f:
        lines = f.readlines()

    # Extract box dimensions from the last line
    box_line = lines[-1].strip().split()
    box = [float(coord) for coord in box_line]

    # Extract coordinates
    coords = []
    atom_names = []
    for line in lines[2:-1]:  # Skip the first two lines and the last line
        parts = line.split()
        coords.append([float(parts[3]), float(parts[4]), float(parts[5])])
        atom_names.append(parts[1])

    return np.array(atom_names), np.array(coords), np.array(box)


def parse_hbond_ndx(ndx_file: str) -> list:
    """Parses a GROMACS index file to extract hydrogen bond indices. The
    index file is expected to contain a section starting with "[ hbonds_"
    and containing lines with three integers representing the donor,
    hydrogen, and acceptor atom indices.
    """
    hbond_idx = list()

    ndx_path = Path(ndx_file)
    if not ndx_path.is_file():
        raise FileNotFoundError(f"File {ndx_file} does not exist.")
    with ndx_path.open("r") as f:
        for line in f:
            if "[ hbonds_" in line:
                break

        for line in f:
            try:
                idx = [int(val) for val in line.strip().split()]
                hbond_idx.append(idx)
            except:
                break

    return hbond_idx


def hbond_data(xpm_file: str, ndx_file: str, gro_file: str):
    """
    Reads hydrogen bond data from an XPM file, index file, and GROMACS .gro file,
    and prints a summary of hydrogen bonds.
    """

    atom_names, _, _ = read_gro(gro_file)
    hbond_idx = parse_hbond_ndx(ndx_file)
    hbond_matrix = parse_xpm_hbond_matrix(xpm_file)

    hbond_sum = hbond_matrix.sum(axis=1)

    nbonds, nframes = hbond_matrix.shape
    if nbonds != len(hbond_idx):
        raise ValueError(
            f"Number of bonds in index file ({len(hbond_idx)}) does not match matrix size ({nbonds})."
        )

    print(f"{'donor':>10} {'H-donor':>10} {'acceptor':>10}")
    print("-" * 33)
    for i, (iat, jat, kat) in enumerate(hbond_idx):
        bond = f"{atom_names[iat-1]:>5}({iat:3d}) {atom_names[jat-1]:>5}({jat:3d}) {atom_names[kat-1]:>5}({kat:3d})"
        line = bond
        line += f" {hbond_sum[i]:5d} {hbond_sum[i] / nframes * 100:>7.2f} %"
        print(line)


if __name__ == "__main__":
    args = sys.argv[1:]
    print(f"Arguments: {args}")

    hbond_data(xpm_file=args[0], ndx_file=args[1], gro_file=args[2])
