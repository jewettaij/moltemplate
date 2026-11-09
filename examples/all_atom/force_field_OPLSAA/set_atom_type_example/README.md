set atom example
====================
This contrived example demonstrates how you can use the "set atom" command to transmute one atom type to another.  When you place a "set atom" inside a "write_once('In Set Atom')" section, bonded force-field parameters for all nearby atoms will be updated.  (Otherwise, they won't.)

This example was taken from the moltemplate manual.


There are no "run" files associated with this example.


### Instructions

More detailed instructions on how to build LAMMPS input files and
run a short simulation are provided in other README files.

To build the files which LAMMPS needs, follow the instructions in:
[README_setup.sh](README_setup.sh)

(The instructions in "README_remove_irrelevant_info.sh" are optional.)







### Customizing atomic charges

In this example, atomic charge for OPLSAA atoms is determined by @atom type
*(...according to a lookup table located at the beginning of the
["oplsaa2024.lt"](../../../moltemplate/force_fields/oplsaa2024.lt) file)*.
*(Any atomic charges listed in the "Data Atoms" section of your molecules'
LT files will be ignored.)*
**These charges can be overridden.**
See [here](../README.md#Customizing-atomic-charges-in-OPLSAA-molecules)
for instructions explaining how to customize atomic charge.
