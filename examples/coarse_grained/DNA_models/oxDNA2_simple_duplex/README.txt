
This is a simple example demonstrating how to build a short straight
coarse-grained helical double-stranded DNA polymer
using the oxDNA2 force field with moltemplate.

This is based on the corresponding LAMMPS example located in this subdirectory:
examples/USER/cgdna/examples/oxDNA2/duplex1/
...which is distributed with LAMMPS.

WARNING: As of 2019-11-02, this example has not been tested or optimized.
         You should probably alter with the timestep, and langevin-settings
         to improve the simulation efficiency.

----
Note: You must compile LAMMPS with the following optional packages installed:
            user-cgdna  asphere
      To do this, go to the "src" directory of your lammps installation and type
      make yes-asphere
      make yes-user-cgdna
      make clean-all
      make NAME_OF_TARGET    # <-- (eg. "make ubuntu", "make g++", "make linux")
----

Instructions on how to build LAMMPS input files and
run a short simulation are provided in other README files.

step 1)
README_setup.sh

step2)
README_run.sh
