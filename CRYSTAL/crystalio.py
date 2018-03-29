#!/usr/bin/env python3
# -*- coding=utf-8 -*-

"""
crystailio
----------

Read a CRYSTAL output file and export all structures, SCF energies and
convergence data.
"""

__author__ = "Germain Vallverdu"
__email__ = "germain.vallverdu@univ-pau.fr"
__licence__ = "GPL"

import re
import numpy as np
from pymatgen import Structure, Lattice


class CrystalOutfile:
    """
    A convenient parser for CRYSTAL output files.

    Args:
        filename: Filename of CRYSTAL output file.
        encoding: encoding of the output file (utf-8)

    Attributes:

    .. attribute:: structures

        All structures from the calculation in the standard orientation. If the
        symmetry is not considered, the standard orientation is not printed out
        and the input orientation is used instead. Check the `standard_orientation`
        attribute.


    """
    def __init__(self, filename, encoding="utf-8"):
        self.filename = filename
        self.encoding = encoding
        self._parse()

    @property
    def initial_structure(self):
        """ First structure read in the calculation """
        return self.structures[0]

    @property
    def final_structure(self):
        """ Last structure read in the calculation """
        return self.structures[-1]

    def get_structure(self, n):
        """ return the nth structure read in the calculation """
        if n > len(self.structures) or n < 1:
            raise ValueError("Bad n value. n=%d. I read %d structures." % (n, len(self.structures)))
        return self.structures[n - 1]

    @property
    def final_energy(self):
        return self.energies[-1]

    def _parse(self):

        float_patt = re.compile(r"[+-]?\d+\.\d+[EFD]?[+-]?\d+")  # -9.3892E+02
        start_patt = re.compile(r"^\s*EEEEEEEEEE STARTING  DATE \d+")
        coord_patt = re.compile(r"^\s+(\d+)\s+(?P<aunit>[TF])\s+(?P<Z>\d+)\s+"
                                r"(?P<specie>\w+)\s+(?P<x>[+-]?\d+\.\d+E[+-]?\d+)"
                                r"\s+(?P<y>[+-]?\d+\.\d+E[+-]?\d+)\s+"
                                r"(?P<z>[+-]?\d+\.\d+E[+-]?\d+)")
        coord_nanotube_patt = re.compile(r"^\s+(\d+)\s+(?P<aunit>[TF])\s+(?P<Z>\d+)\s+"
                                         r"(?P<specie>\w+)\s+(?P<x>[+-]?\d+\.\d+E[+-]?\d+)"
                                         r"\s+(?P<y>[+-]?\d+\.\d+E[+-]?\d+)\s+"
                                         r"(?P<z>[+-]?\d+\.\d+E[+-]?\d+)\s+"
                                         r"(?P<radius>\d+\.\d+)")
        forces_patt = re.compile(r"^\s+(?P<iat>\d+)\s+(?P<Z>\d+)\s+"
                                 r"(?P<x>[+-]?\d+\.\d+E[+-]?\d+)\s+"
                                 r"(?P<y>[+-]?\d+\.\d+E[+-]?\d+)\s+"
                                 r"(?P<z>[+-]?\d+\.\d+E[+-]?\d+)")
        max_grad_patt = re.compile(r"^\sMAX GRADIENT\s+(?P<max_grad>\d+\.\d+)"
                                   r"\s+THRESHOLD\s+(?P<max_grad_thr>\d+\.\d+)")
        rms_grad_patt = re.compile(r"^\sRMS GRADIENT\s+(?P<rms_grad>\d+\.\d+)"
                                   r"\s+THRESHOLD\s+(?P<rms_grad_thr>\d+\.\d+)")
        max_displac_patt = re.compile(r"^\sMAX DISPLAC\.\s+(?P<max_displac>\d+\.\d+)"
                                      r"\s+THRESHOLD\s+(?P<max_displac_thr>\d+\.\d+)")
        rms_displac_patt = re.compile(r"^\sRMS DISPLAC\.\s+(?P<rms_displac>\d+\.\d+)"
                                      r"\s+THRESHOLD\s+(?P<rms_displac_thr>\d+\.\d+)")
        norm_grad_patt = re.compile(r"^\s+GRADIENT NORM\s+(?P<norm_grad>\d+\.\d+)"
                                    r"\s+GRADIENT THRESHOLD\s+(?P<norm_grad_thr>\d+\.\d+)")

        self.title = ""
        self.system = ""
        self.group = ""
        self.slab = False
        self.nanotube = False
        self.volumes = list()
        self.energies = list()
        self.forces = list()
        self.convergence_data = list()
        self.geometry_converge = False
        self.scf_converge = False

        external_geometry = False
        with open(self.filename, "r", encoding=self.encoding) as f:

            # look for starting message
            for line in f:
                if start_patt.match(line):
                    self.title = f.readline().strip()
                    break

            # ------------------------------------------------------------------
            # first, read the initial geometry & identify the type of structure
            # ------------------------------------------------------------------
            for line in f:
                if re.match(r"^\sGEOMETRY INPUT FROM EXTERNAL FILE", line):
                    external_geometry = True
                    line = f.readline()
                    if "SLAB" in line:
                        self.slab = True
                    if "NANOTUBE" in line:
                        self.nanotube = True

                    print("WARNING: Geometry from an external file.")
                    break

                if re.match(r"^\sSLAB CALCULATION", line):
                    self.slab = True
                    system_patt = re.compile(r"^\sSYSTEM AND LATTICE")
                    group_patt = re.compile(r" PLANE GROUP N.")
                    break

                if re.match(r"^\sCRYSTAL CALCULATION", line):
                    system_patt = re.compile(r"^\sCRYSTAL FAMILY")
                    group_patt = re.compile(r"^\sSPACE GROUP")
                    break

            # look for initial geometry: GEOMETRY FOR WAVEFUNCTION
            # read group and crystallographic system
            # check if a SLAB or NANOTUBE is built by GEOMETRY EDITING
            geom_for_wf = False
            for line in f:
                if not external_geometry and system_patt.search(line):
                    self.system = line.split(":")[1].strip()

                if not external_geometry and group_patt.search(line):
                    self.group = line.split(":")[1].strip()

                if " SLAB GENERATED " in line:
                    self.slab = True
                    # group and system no more relevant
                    self.group = ""
                    self.system = ""

                if "CONSTRUCTION OF A NANOTUBE FROM A SLAB" in line:
                    self.nanotube = True
                    self.slab = False
                    # group and system no more relevant
                    self.group = ""
                    self.system = ""

                if re.match(r"^\sGEOMETRY FOR WAVE FUNCTION", line):
                    geom_for_wf = True
                    break

                if line == "":
                    # end of file, geometry for wavefunction not found
                    break

            if not geom_for_wf:
                # STOP case, add TESTGEOM to d12
                raise ValueError("GEOMETRY FOR WAVEFUNCTION NOT FOUND.\n"
                                 "Please, add TESTGEOM in the d12 input file.")

            # read until calculation start
            # read starting geometry and look for PRIMITIVE or CRYSTALLOGRAPHIC
            read_geom = False
            while "CRYSTAL - SCF - TYPE OF CALCULATION" not in line:
                line = f.readline()

                if line == "":
                    raise ValueError("End of file.")

                # search PRIMITIVE CELL
                if re.match(r"^\sPRIMITIVE CELL", line):
                    read_geom = True
                    geom_patt = re.compile(r"^\sPRIMITIVE CELL")

                # search CRYSTALLOGRAPHIC CELL if exist
                if re.match(r"^\sCRYSTALLOGRAPHIC CELL", line):
                    read_geom = True
                    geom_patt = re.compile(r"^\sCRYSTALLOGRAPHIC CELL")

                if read_geom:
                    if not self.slab and not self.nanotube:
                        volume = float(line.split("=")[1].split()[0].strip(")"))
                        self.volumes.append(volume)
                    f.readline()

                    # lattice parameters
                    line = f.readline()
                    params = [float(val) for val in re.findall(r"\d+\.\d+", line)]
                    lattice = Lattice.from_lengths_and_angles(params[0:3], params[3:])

                    # step on for 4 lines
                    [f.readline() for _ in range(4)]

                    # read coordinates
                    species = list()    # atom names
                    uniq = list()   # True if atom belong to the asymmetric unit
                    radius = list()  # distance from the axes of the nanotube
                    coords = list()
                    while line != "\n":
                        read = False
                        line = f.readline()
                        if self.nanotube and coord_nanotube_patt.match(line):
                            data = coord_nanotube_patt.match(line).groupdict()
                            read = True
                        elif coord_patt.match(line):
                            data = coord_patt.match(line).groupdict()
                            read = True

                        if read:
                            specie = data["specie"]
                            specie = specie if len(specie) == 1 else specie[0] + specie[1].lower()
                            species.append(specie)
                            coord = [float(data[k]) for k in "xyz"]
                            uniq.append(True if data["aunit"] == "T" else False)

                            if self.slab:
                                coord[2] /= lattice.c
                            elif self.nanotube:
                                coord[1] /= lattice.b
                                coord[2] /= lattice.c
                                radius.append(float(data["radius"]))

                            coords.append(coord)

                    self.structures = [Structure(lattice, species, coords,
                                                 site_properties={"aunit": uniq})]

                    read_geom = False

            # ------------------------------------------------------------------
            # from that point, SCF, or structure optimization start !
            # continue up to the end of file
            # ------------------------------------------------------------------
            n_geom = 0
            cvg_data = dict()
            while line != "":
                line = f.readline()

                if " TOTAL ENERGY" in line:
                    self.energies.append(float(float_patt.findall(line)[0]))
                    self.scf_converge = True

                if "CARTESIAN FORCES IN HARTREE/BOHR" in line:
                    # WARNING: Forces are not printed at each geom step
                    line = f.readline()
                    forces = list()
                    for _ in range(self.initial_structure.num_sites):
                        data = forces_patt.match(f.readline()).groupdict()
                        forces.append([float(data[c]) for c in "xyz"])
                    self.forces.append(np.array(forces))

                if max_grad_patt.match(line):
                    cvg_data.update(max_grad_patt.match(line).groupdict())
                if rms_grad_patt.match(line):
                    cvg_data.update(rms_grad_patt.match(line).groupdict())
                if max_displac_patt.match(line):
                    cvg_data.update(max_displac_patt.match(line).groupdict())
                if rms_displac_patt.match(line):
                    cvg_data.update(rms_displac_patt.match(line).groupdict())
                if norm_grad_patt.match(line):
                    cvg_data.update(norm_grad_patt.match(line).groupdict())

                if line == "":
                    # end of file ?
                    break

                if "COORDINATE AND CELL OPTIMIZATION" in line:
                    cvg_data = {k: float(v) for k, v in cvg_data.items()}
                    # end of optimization cycle
                    self.convergence_data.append(cvg_data)
                    n_geom += 1
                    cvg_data = dict()

                if "FINAL OPTIMIZED GEOMETRY" in line:
                    self.geometry_converge = True
                    n_geom += 1

                # search structure data
                if geom_patt.match(line):
                    # PRIMITVE or CRYSTALLOGRAPHIC depending on what is present
                    read_geom = True

                if read_geom:
                    if not self.slab and not self.nanotube:
                        volume = float(line.split("=")[1].split()[0].strip(")"))
                        self.volumes.append(volume)
                    f.readline()

                    # lattice parameters
                    line = f.readline()
                    params = [float(val) for val in re.findall(r"\d+\.\d+", line)]
                    lattice = Lattice.from_lengths_and_angles(params[0:3], params[3:])

                    # step on for 4 lines
                    [f.readline() for _ in range(4)]

                    # read coordinates
                    species = list()    # atom names
                    uniq = list()   # True if atom belong to the asymmetric unit
                    radius = list()  # distance from the axes of the nanotube
                    coords = list()
                    while line != "\n":
                        read = False
                        line = f.readline()
                        if self.nanotube and coord_nanotube_patt.match(line):
                            data = coord_nanotube_patt.match(line).groupdict()
                            read = True
                        elif coord_patt.match(line):
                            data = coord_patt.match(line).groupdict()
                            read = True

                        if read:
                            specie = data["specie"]
                            specie = specie if len(specie) == 1 else specie[0] + specie[1].lower()
                            species.append(specie)
                            coord = [float(data[k]) for k in "xyz"]
                            uniq.append(True if data["aunit"] == "T" else False)

                            if self.slab:
                                coord[2] /= lattice.c
                            elif self.nanotube:
                                coord[1] /= lattice.b
                                coord[2] /= lattice.c
                                radius.append(float(data["radius"]))

                            coords.append(coord)

                    self.structures.append(Structure(lattice, species, coords,
                                                     site_properties={"aunit": uniq}))

                    read_geom = False


if __name__ == "__main__":
    pass
