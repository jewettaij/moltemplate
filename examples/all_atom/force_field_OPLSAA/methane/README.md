Methane Example
===================
This example demonstrates how to build a simulation containing a box of methane gas.  (Not a very interesting example.)

### Details 

The methane molecules in this example use the OPLSAA force-field.
This means that the database of force-field parameters in "oplsaa.lt"
will be used to generate angles, dihedrals, and impropers.
The "moltemplate_files/methane.lt" file
contains these lines which refer to OPLSAA:

```
import "oplsaa.lt"
Methane inherits OPLSAA { ... }    # (see "methane.lt")
```


### Instructions

More detailed instructions on how to build LAMMPS input files and
run a short simulation are provided in other README files.

1) To build the files which LAMMPS needs, follow the instructions in:
[README_setup.sh](README_setup.sh)

2) To run LAMMPS with these files, follow these instructions:
[README_run.sh](README_run.sh)

*(If you notice a problem with this example, please [report it](../README.md).)*
