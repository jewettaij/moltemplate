Benzene Example
===================
This simple example demonstrates how to build a simulation containing a box full of benzene molecules and run a simulation at conditions of 298K, 1 barr.  The number of benzene molecules, and box size can be controlled by editing the [system.lt](./moltemplate_files/system.lt) file.  The simulation conditions can be controlled by editing the [run.in.npt](./run.in.npt) file.


#### Images

<img src="images/benzene.jpg" width=110>

The number of molecules in the simulation and and the simulation box size can be controlled by editing the [system.lt file](moltemplate_files/system.lt).  The simulation contitions can be controlled by editing the [run.in.npt file](run.in.npt).


### Details 

The benzene molecules in this example use the AMBER (GAFF2) force-field.  *(The GAFF2 force-field is also available.)*  This means that the database of force-field parameters in "gaff2.lt" will be used to generate angles, dihedrals, and impropers.  The "moltemplate_files/benzene.lt" file contains these lines which refer to GAFF2:

```
import "gaff2.lt"
Benzene inherits GAFF2 { ... }    # (see "benzene.lt")
```


### *WARNING: The atomic charges in this examples are not correct*

The AMBER for field does not include charge information.  (In this example, they were borrowed from the corresponding atoms in the ["oplsaa2023.lt" file](../../../../moltemplate/force_fields/oplsaa2023.lt).  Do not do this!)

The generation of atomic partial charges requires 3rd party software.

For suggestions how to calculate charges correctly, see [README.md](../README.md).  Then choose a 3rd party program to calculate partial charges for the atoms in each type of molecule.  Then edit the corresponding ".lt" files in the "moltemplate_files" directory accordingly.


### Instructions

More detailed instructions on how to build LAMMPS input files and
run a short simulation are provided in other README files.

1) To build the files which LAMMPS needs, follow the instructions in:
[README_setup.sh](README_setup.sh)

2) To run LAMMPS with these files, follow these instructions:
[README_run.sh](README_run.sh)

(The instructions in "README_remove_irrelevant_info.sh" are optional.  *(If you notice a problem with this example, please [report it](../README.md).*)


### Requirements

This example requires a version of LAMMPS compiled with support for the optional "EXTRA-MOLECULE" package (because the AMBER force field currently uses *dihedral_style fourier*).  If you encounter the error *"Invalid dihedral_style"*, then see [this page](https://docs.lammps.org/Build_package.html) for instructions to compile LAMMPS to support this package.
