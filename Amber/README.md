# AMBER : input files, scripts, fortran topology module

This folder contains **OLD** [Amber](http://ambermd.org/) input files
or scripts (2010).

## Read topology file

This is a Fortran90 module which contains variables for

* The atom number and the residu number
* A vector containing the residu number of each atom
* The number of the first atom of each residu
* The mass and the charge of each atom
* The name of each atom
* The residu name of each atom
* The list of residu


There is also a subroutine which allocate each vector or table and read the topology file.

## AMBER input files and Scripts

Here are input files to run classical molecular dynamics or minimization with sander or pmemd program of [Amber](http://ambermd.org/).
Scripts files are bash commands which generate the MDIN file and other necessary input file and run the calculation.

MDIN file are classical sander/pmemd input files. Each line starting with \# are
commentary and must be removed before you used it.

The following MDIN files are given in a specific order.
You can follow this order to run (N,P,T) simulation from
a pdb file of the [Protein Data Bank](http://www.rcsb.org/)
(you have to prepare topology and coordinate files before you start the first minimization).

* **minH.in**: Minimization of all hydrogen atoms. Constraints are put on all heavy atoms.
* **min_wat.in**: Minimization of all water molecules. Constraint are put on the protein.
* **mintot.in**: Total minimization of the protein and water molecules with small constraints on the protein.
* **eq-NVT.in**: Simulation at constant volume with small constraints on the protein. Now we warm slowly the simulation box from 100K to 300K.
* **sim-NPT.in**: (N,P,T) Simulation at T = 300 K and P = 1.0 atm.
* **sim-NVE.in**: (N,V,E) Simulation.
* **bias_MD_XX.in**: Bias simulations. You also nead forget the **disang** file in
order to define the bias potential.
