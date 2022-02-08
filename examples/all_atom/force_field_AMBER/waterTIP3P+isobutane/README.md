Isobutane and water phase separation example
==============
A mixture of two small organic molecules using the *AMBER/GAFF* force field.  In this example, the water molecules were initially arranged in a rectangular lattice.  The isobutane molecules were also arranged in a lattice, and were shifted to avoid overlap with the water molecules. *(Alternatively, you can create a single lattice and specify the number of isobutane and water molecules you want in it using moltemplate's "new random([],[])" command, which is explained in [the manual](https://moltemplate.org/doc/moltemplate_manual.pdf#subsubsection.8.9.1).  This gives you more control over the concentration of each ingredient.  You can also use PACKMOL to create random mixtures of molecules.)*  The two types of molecules phase separate over the course of a few hundred ps.  The GAFF parameters are applied only to the isobutane molecule.  (The water molecule paramters are defined explicitly in the "force_fields/tip3p_2004.lt" file distributed with moltemplate.)


#### Images

<img src="images/isobutane.jpg" width=110> <img src="images/plus.svg" height=80> <img src="images/water.jpg" width=110> <img src="images/rightarrow.svg" height=80> <img src="images/water+isobutane_t=0_LR.jpg" width=150> <img src="images/rightarrow.svg" height=80> <img src="images/water+isobutane_t=840ps_LR.jpg" width=150>

The number of molecules, positions, and simulation box size can be controlled by editing the [system.lt file](moltemplate_files/system.lt).  The simulation contitions can be controlled by editing the [run.in.npt file](run.in.npt).


### *WARNING: The atomic charges in this examples are not correct*

The AMBER for field does not include charge information.  (In this example, they were borrowed from the corresponding atoms in the ["oplsaa.lt" file](../../../../moltemplate/force_fields/oplsaa.lt).  Do not do this!)

The generation of atomic partial charges requires 3rd party software.

For suggestions how to calculate charges correctly, see [README.md](../README.md).  Then choose a 3rd party program to calculate partial charges for the atoms in each type of molecule.  Then edit the corresponding ".lt" files in the "moltemplate_files" directory accordingly.


## Instructions

1) To build the files which LAMMPS needs, follow the instructions in:
[README_setup.sh](README_setup.sh)

2) To run LAMMPS with these files, follow these instructions:
[README_run.sh](README_run.sh)

(The instructions in "README_remove_irrelevant_info.sh" are optional.  *(If you notice a problem with this example, please [report it](../README.md).*)


### Requirements

This example requires a version of LAMMPS compiled with support for the optional "USER-MISC" package (because the AMBER force field currently uses *dihedral_style fourier*).  If you encounter the error *"Invalid dihedral_style"*, then see [this page](https://lammps.sandia.gov/doc/Build_package.html) for instructions to compile LAMMPS with the "USER-MISC" package.
