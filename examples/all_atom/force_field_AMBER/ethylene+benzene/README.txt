This example shows how to simulate a mixture of ethylene and benzene
using the AMBER/GAFF force field.

WARNING:
The atomic partial charges in this example are not correct!
For details how to calculate charges correctly, see:
https://github.com/jewettaij/moltemplate/blob/master/examples/all_atom/force_field_AMBER/README.md
...and edit the .lt files in the moltemplate_files directory accordingly.

As of 2016-11-21, this code has not been tested for accuracy.
(See the WARNING.TXT file.)

step 1)
To build the files which LAMMPS needs, follow the instructions in:
README_setup.sh

step 2)
To run LAMMPS with these files, follow these instructions:
README_run.sh
