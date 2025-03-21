Hexadecane example
==============
This example is a simple simulation of many long alkane chains (hexadecane) in a box near the boiling point at atmospheric pressure.  The hexadecane molecule in this example (defined in the [hexadecane.lt](moltemplate_files/hexadecane.lt) file) was constructed from monomeric subunits (named "CH2", and "CH3").

#### Images

<img src="images/ch2_ry60_LR.jpg" width=110> <img src="images/plus.svg" height=80> <img src="images/ch3_ry60_LR.jpg" width=110> <img src="images/rightarrow.svg" height=80> <img src="images/hexadecane_LR.jpg" width=150>  <img src="images/rightarrow.svg" height=80> <img src="images/hexadecane_12x12x2_t=0_LR.jpg" width=150> <img src="images/rightarrow.svg" height=80> <img src="images/hexadecane_12x12x2_t=10ps_npt_LR.jpg" width=150>

The number of molecules and simulation box size can be controlled by editing the [system.lt file](moltemplate_files/system.lt).  The length of each polymer can be controlled by editing the [hexadecane.lt](moltemplate_files/hexadecane.lt) file.  The simulation contitions can be controlled by editing the [run.in.npt file](run.in.npt).


### *WARNING: The atomic charges in this examples are not correct*

The AMBER for field does not include charge information.  (In this example, they were borrowed from the corresponding atoms in the ["oplsaa2024.lt" file](../../../../moltemplate/force_fields/oplsaa2024.lt).  Do not do this!)

The generation of atomic partial charges requires 3rd party software.

For suggestions how to calculate charges correctly, see [README.md](../README.md).  Then choose a 3rd party program to calculate partial charges for the atoms in each type of molecule.  Then edit the corresponding ".lt" files in the "moltemplate_files" directory accordingly.


## Instructions

1) To build the files which LAMMPS needs, follow the instructions in:
[README_setup.sh](README_setup.sh)

2) To run LAMMPS with these files, follow these instructions:
[README_run.sh](README_run.sh)

(The instructions in "README_remove_irrelevant_info.sh" are optional.  *(If you notice a problem with this example, please [report it](../README.md).*)


### Requirements

This example requires a version of LAMMPS compiled with support for the optional "EXTRA-MOLECULE" package (because the AMBER force field currently uses *dihedral_style fourier*).  If you encounter the error *"Invalid dihedral_style"*, then see [this page](https://docs.lammps.org/Build_package.html) for instructions to compile LAMMPS to support this package.
