Examples
================
This directory contains examples for the
"moltemplate" molecule builder for LAMMPS.
http://www.moltemplate.org

Each directory contains one or more examples.

Each example directory contains the following files and directories:

| File or directory   | Explanation
|---------------------|--------------------------------------------------------|
| images/             | This folder has pictures of the molecules in the system|
| moltemplate_files/  | This folder contains LT files and other auxiliary files|
| README_setup.sh     | Instructions for how to use moltemplate to create files|
| README_run.sh       | Instructions for how to run LAMMPS (after moltemplate) |
| README_visualize.txt| Instructions for viewing in DATA/DUMP files in VMD     |

...and one or more LAMMPS input scripts with names like "run.in.min",
"run.in.npt", and "run.in.nvt". These input scripts load the files you
created with moltemplate which describe the system that you wish to simulate.
They also contain a (usually minimal) list of additional LAMMPS
commands needed to run the simulation under reasonable conditions
(specifying timesteps, integrators, and output file names, for example).

You can run LAMMPS using commands like this:
```
   lmp_mpi -i run.in.npt
```
(See the "README_run.sh" file in each directory for details.
 The name of your lammps binary, "lmp_mpi" in this example, may vary.
 Sometimes, these scripts must be run in a certain order.  For example
 it may be necessary to run "run.in.min" to minimize the system before
 you can use "run.in.npt", and later "run.in.nvt".  Many of these
 script files have not been optimized.)

