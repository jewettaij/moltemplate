Methane Example
===================
This example demonstrates how to build a simulation containing a box of methane gas.  (Not a very interesting example.)


### Details 

The methane molecules in this example use the OPLSAA force-field.  This means that the database of force-field parameters in "oplsaa.lt" will be used to generate angles, dihedrals, and impropers.  The "moltemplate_files/methane.lt" file contains these lines which refer to OPLSAA:

```
import "oplsaa.lt"
Methane inherits OPLSAA { ... }    # (see "methane.lt")
```

The number of molecules and simulation box size can be controlled by editing the [system.lt file](moltemplate_files/system.lt).  The simulation contitions can be controlled by editing the [run.in.npt file](run.in.npt).


### Instructions

More detailed instructions on how to build LAMMPS input files and run a short simulation are provided in other README files.

1) To build the files which LAMMPS needs, follow the instructions in:
[README_setup.sh](README_setup.sh)

2) To run LAMMPS with these files, follow these instructions:
[README_run.sh](README_run.sh)

(The instructions in "README_remove_irrelevant_info.sh" are optional.  *(If you notice a problem with this example, please [report it](../README.md).*)


### Customizing atomic charges

In this example, atomic charge for OPLSAA atoms is determined by @atom type
*(...according to a lookup table located at the beginning of the
["oplsaa.lt"](../../../moltemplate/force_fields/oplsaa.lt) file)*.
*(Any atomic charges listed in the "Data Atoms" section of your molecules'
LT files will be ignored.)*
**These charges can be overridden.**
See [here](../README.md#Customizing-atomic-charges-in-OPLSAA-molecules)
for instructions explaining how to customize atomic charge.
