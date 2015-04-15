#!/usr/bin/env python
# -*- coding=utf-8 -*-

from __future__ import print_function, division

"""
doc
"""
__author__ = "Germain Vallverdu <germain.vallverdu@univ-pau.fr>"
__licence__ = "GPL"
__date__ = "Février 2015"

import numpy as np
import matplotlib.pyplot as plt
import pymatgen as mg

#-----------------------------------------------------------------------------
# ideal molecule constructor
#-----------------------------------------------------------------------------
def mol_Oh(central, ligand, scale):
    """ 
    Return a perfect octahedra as a mg.Molecule object.

    Args:
        central: (string) Name of the central atom
        ligand: (string) Name of ligand atoms
        scale: (float) length of central-ligand distance

    Returns
        mg.Molecule object
    """
    species = [mg.Element(central)] + 6 * [mg.Element(ligand)]
    template = [[ 0., 0., 0.],
                [ 1., 0., 0.],
                [-1., 0., 0.],
                [ 0., 1., 0.],
                [ 0.,-1., 0.],
                [ 0., 0., 1.],
                [ 0., 0.,-1.]]
    coords = [[scale * xi for xi in coord] for coord in template]
    return mg.Molecule(species, coords)

def mol_D4h(central, ligand, d1, d2):
    """ 
    Return a D4h molecule as a mg.Molecule object.

    Args:
        central: (string) Name of the central atom
        ligand: (string) Name of ligand atoms
        d1: (float) length of central-ligand axial distance
        d2: (float) length of central-ligand equatorial distance

    Returns
        mg.Molecule object
    """
    species = [mg.Element(central)] + 6 * [mg.Element(ligand)]
    template = [[ 0., 0., 0.],
                [ 1., 0., 0.],
                [-1., 0., 0.],
                [ 0., 1., 0.],
                [ 0.,-1., 0.],
                [ 0., 0., 1.],
                [ 0., 0.,-1.]]
    coords = [template[0]] \
           + [[d1 * xi for xi in coord] for coord in template[1:5]] \
           + [[d2 * xi for xi in coord] for coord in template[5:]]
    return mg.Molecule(species, coords)

#-----------------------------------------------------------------------------
# Distortion measurment
#-----------------------------------------------------------------------------
def centralDist(envs, ideal, dmin=1., dmax=3.5, step=0.05):
    """ 
    Plot an histogram of central to ligand distances between dmin and dmax with
    a bin size step.

    Args:
        envs: list of mg.Molecule object. The first atom is the central atom
        the following atoms are the ligand.
        dmin (float): lower bound of the histogram
        dmax (float): upper bound of the histogram
        step (float): bin size

    """ 
    # compute all central - ligand distances
    data = list()
    for mol in envs:
        data += sorted([di for n, di in mol.get_neighbors(mol[0], 10.)])

    # print data
    nligand = len(envs[0][1:])
    nsite = len(envs)
    np.savetxt("all_env.dat", 
               np.array(data).reshape(nsite, nligand).transpose(), 
               fmt="%10.4f",
               header="Each column contains the sorted central-ligand" \
               + " distances of one site")

    # plot data
    font = {'family': 'serif', 'size': 20}
    plt.rc('font', **font)                                                           
    plt.figure(figsize=(11.69, 8.27)) # A4

    mybins = np.arange(dmin, dmax, step)
    hval, bins, p = plt.hist(data, mybins, normed=False)
    hval /= nsite

    ideal_d = [di for n, di in ideal.get_neighbors(ideal[0], 10.)]
    hid, b, p = plt.hist(ideal_d, mybins, normed=False)

    plt.clf()

    plt.bar(bins[:-1], hval, width=step, color="blue", label="structure")
    plt.bar(bins[:-1], hid, width=step, color="green", label="reference")
    plt.xlabel("distances (A)")
    plt.ylabel("Histogram")
    plt.legend()
    plt.grid()
    plt.xlim(dmin, dmax)
    plt.show()

def paramVector(mol, radius = 10.):
    """ 
    Build the vector representing the site in the parameters' space.

    Args:
        * mol : mg.Molecule object of the central atom plus the ligand atoms.
              : The first atom must be the central atom
        * radius : cut off distance for neighbors of the central atoms

    Returns:
        * a numpy array. The n first distances are the central to all ligand
        distances sorted from the smallest to the largest. The m following
        distances are all ligand - ligand distances sorted from the smallest to
        the largest.
    """
    # distances central atom to ligand
    d_central = [di for n, di in mol.get_neighbors(mol[0], radius)]
    # all distances
    dm = mol.distance_matrix
    d_all = [di for line in np.triu(dm) for di in line if di > 0.]
    return np.array(sorted(d_central) + sorted(d_all))

def readStructures(sfile, central, ligand, radius=3.):
    """
    Read structure from file sfile and look for environments of all central
    atoms in the structure. Only neighbors of type ligand are taken into
    account.

    Returns a list of mg.Molecule in which the first atom is the central atom
    and the following atoms are the neighbors of type ligand of this central
    atom.
    """

    struct = mg.Structure.from_file(sfile)
    central_sites = [site for site in struct 
                          if site.specie == mg.Element(central)]

    envs = list()

    print("Identified sites:")
    for site in central_sites:

        mol = mg.Molecule([site.specie], [site.coords])

        neighbors = struct.get_neighbors(site, radius)

        for neighbor, d in neighbors:
            if neighbor.specie == mg.Element(ligand):
                mol.append(neighbor.specie, neighbor.coords)

        envs.append(mol)
        print("%2s[%2d] : #ligand = %d (%s)" %
              (site.specie, 
               struct.index(site), 
               len(mol[1:]), 
               " ".join([at.specie.symbol for at in mol[1:]])))

    return envs

def norm(vec, origin):
    """ 
    return the norm considering the vector vec in the parameters' space and
    the ideal vector origin in the parameters' space.
    """
    return np.sqrt(((vec - origin)**2).sum()/(origin**2).sum())

if __name__ == "__main__":

    # ideal
    #envfile = "ideal.xyz"
    #ideal_mol = mg.Molecule.from_file(envfile)
    #ideal_mol = mol_Oh("Mn", "O", 2.)
    ideal_mol = mol_D4h("Mn", "O", 1.9, 2.2)
    ideal_mol.to("xyz", "d4h.xyz")
    idealVec = paramVector(ideal_mol)

    # analysed structure
    radius = 2.5
    sfile = "POSCAR"
    envs = readStructures(sfile, "Mn", "O", radius)

    for mol in envs:
        vec = paramVector(mol)
        print(norm(vec, idealVec))

    centralDist(envs, ideal_mol)

