#!/usr/bin/env python
# -*- coding=utf-8 -*-

from pymatgen.core import Molecule, SymmOp
import numpy as np


def local_frame(mol, origin=0, xat=1, yat=2, barycenter=False):
    """
    Return a new molecule object where the atoms are in the local frame defined
    from the three atoms given as argument.

    Args:
        mol (Molecule): The molecule to adsorb.
        origin (int): origin of the frame
        xat (int): atom at the end of x axis
        yat (int): third atom to define the frame

    Retuns
        Returns a new molecule object with coordinates in the frame.

    """

    # set up the local frame
    u = mol[xat - 1].coords - mol[origin - 1].coords
    u /= np.linalg.norm(u)
    v = mol[yat - 1].coords - mol[origin - 1].coords
    v = v - np.dot(u, v) * u
    v /= np.linalg.norm(v)
    w = np.cross(u, v)

    mat = np.array([u, v, w])

    # move molecule around an origin
    G = np.zeros(3)
    if barycenter is True:
        G = mol.cart_coords.sum(axis=0) / len(mol)
    elif len(barycenter) == 3:
        G = np.array(barycenter)

    # new coordinates in the local frame
    newcoords = np.dot(mol.cart_coords - G, mat.T)

    return Molecule(mol.species, newcoords, mol.charge, mol.spin_multiplicity,
                    mol.site_properties)


def adsorb(mol, slab, atmol, atslab, z, img=None, shift=[0., 0., 0.]):
    """
    Add molecule `mol` on the surface `slab` at a distance from the slab `z`. The
    list `atmol` defines the atoms of the molecule which will be adsorbed. The list
    `atslab` define the adsorption site on the slab. The function assumes that the slab
    surface is perpendicular to the z axes.

    Args:
        mol (Molecule): The molecule to adsorb.
        slab (Slab): The slab on which adsorption is done.
        atmol ([int]): Sequence of atom's numbers of the molecule which define
            the anchor point of the molecule.
        atslab ([int]): Sequance of atom's numbers of the slab which define the
            adsorption site.
        z (float): distance between the barycenters of atslab and atmol atoms.
        img (Nx3 array): translation vectors which define periodic images of
            atoms belonging to the slab.
        shift ([x, y, z]): add a translation vector to the molecule

    All parameters are mandatory, the periodic images definition (img), the shift
    vector and the verbosity control.

    The adsorption site type on the slab is defined by the `atslab` argument..
    This list contains the atom number of atoms belonging to the slab on which
    we will adsorb the molecule.

        * 1 atom : top site
        * 2 atom : bridge site (2fold site)
        * 3 atom : hollow bridge site (3fold site)
        * 4 atom : 4 fold site
        * ...

    The img argument define for each atom present in the atslab list, the
    periodic image which has to be used in order to build the adsorption site.

    The z axes will be kept perpendicular to the slab. Thus take care to prepare the
    molecule in a way that the adsorption will be done along the z axes. The distance
    z, will be the distance between the atoms of the molecule defined by the parameter
    `atmol` and the barycenter of the slab's atoms defined by the `atslab` parameter. The
    barycenter is not computed using atoms mass.

    """

    # check
    # -----
    if not isinstance(mol, Molecule):
        print("mol must be a Molecule object")
        exit(1)
    if not isinstance(atmol, list):
        print("atmol must be a list of atom number")
        print("atmol : {0}".format(atmol))
        exit(1)
    if img is not None:
        if len(img) != len(atslab):
            print("You must give imaging information for all atslab atoms")
            exit(1)
    if not isinstance(atslab, list):
        print("atslab must be a list of atom number")
        print("atslab : {0}".format(atslab))
        exit(1)
    if not isinstance(shift, list) and not isinstance(shift, np.ndarray):
        print("shift must be a list or a ndarray")
        exit(1)
    elif len(shift) != 3:
        print("shift length must be 3")
        print("len(shift) = " + str(len(shift)))
        exit(1)

    # default value for periodic images
    if img is None:
        img = [[0., 0., 0.] for i in range(len(atslab))]

    # convert shift to numpy array
    if not isinstance(shift, np.ndarray):
        shift = np.array(shift)

    # atom numbers start at 0
    atslab = [iat - 1 for iat in atslab]
    atmol = [iat - 1 for iat in atmol]

    # move slab's atoms of the adsorption site according to img vectors
    # xyz_site are the coordinate of slab's atom which define the asorption site
    xyz_site = list()
    for i, iat in enumerate(atslab):
        imgtrans = np.zeros(3)
        for k in range(3):
            imgtrans += img[i][k] * slab.lattice.matrix[k]
        xyz_site.append(slab[iat].coords + imgtrans)
    xyz_site = np.array(xyz_site)

    # coordinates of the barycenter of the adsorption site : Gslab
    Gslab = xyz_site.mean(axis=0)

    # coordinate of the barycenter of molecule's atoms : Gmol
    Gmol = np.zeros(3)
    for iat in atmol:
        Gmol += mol[iat].coords
    Gmol /= float(len(atmol))

    # vertical vector along which the adsorption is done
    if len(atslab) == 1:
        # top site: adsorption along z axis
        ztrans = slab.normal

    elif len(atslab) == 2:
        # bridge or 2-fold site adsorption perpendicularly to the vector
        # joining the 2 slab atoms
        u = xyz_site[1] - xyz_site[0]
        u /= np.sqrt((u**2).sum())

        ztrans = slab.normal - np.dot(u, slab.normal) * u

    elif len(atslab) == 3:
        # 3-fold site: ztrans is obtained from a cross product
        u = xyz_site[1] - xyz_site[0]
        u /= np.sqrt((u**2).sum())

        v = xyz_site[2] - xyz_site[0]
        v = v - np.dot(u, v) * v
        v /= np.sqrt((v**2).sum())

        ztrans = np.cross(u, v)

        # ztrans and z axes in the same direction
        if np.dot(ztrans, np.array([0., 0., 1.])) < 0:
            ztrans *= -1.

    else:
        # TODO: average plane
        print("WARNING : z shift along z axes from barycenter of atoms")
        ztrans = slab.normal

    # translate the molecule above the adsorption site
    vtrans = Gslab - Gmol + z * ztrans + shift
    trans = SymmOp.from_axis_angle_and_translation(axis=np.array([1., 0., 0.]),
                                                   angle=0.,
                                                   translation_vec=vtrans)
    mol.apply_operation(trans)

    # create the new slab with the adsorbed molecule
    for specie, coords in zip(mol.species, mol.cart_coords):
        slab.append(specie, coords, coords_are_cartesian=True)

    return slab


if __name__ == "__main__":
    import pymatgen as mg
    from pymatgen.core.surface import Slab

    slabfile = "./CONTCAR_1x1.vasp"
    # slabfile = "./CONTCAR"
    s = mg.Structure.from_file(slabfile)
    slab = Slab(s.lattice, s.species, s.frac_coords,
                miller_index=[0, 0, 1],
                oriented_unit_cell=s,
                shift=0.,
                scale_factor=[])

    # ------------------------------
    molxyz = "./h2O.xyz"
    mol = mg.Molecule.from_file(molxyz)
    G = [0, 0, 0]
    mol = local_frame(mol, origin=1, xat=0, yat=2, barycenter=G)

    op = SymmOp.from_origin_axis_angle(mol[1].coords, axis=[0, 1, 0], angle=-90)
    mol.apply_operation(op)

    ads_slab = adsorb(mol, slab.copy(), atmol=[1], atslab=[22], z=0.)
    mg.io.vasp.inputs.Poscar(ads_slab).write_file("ads_H2O.vasp")


    exit()
    # -------------------------------
    # monolayer H2O
    ads_slab = adsorb(mol, slab.copy(), atmol=[2], atslab=[6], z=1.)

    op = SymmOp.from_origin_axis_angle(mol[1].coords, axis=[0, 0, 1], angle=180)
    mol.apply_operation(op)
    ads_slab = adsorb(mol, ads_slab, atmol=[2], atslab=[15], z=-1.)

    mg.io.vasp.inputs.Poscar(ads_slab).write_file("ads1x1_H2O.vasp")


