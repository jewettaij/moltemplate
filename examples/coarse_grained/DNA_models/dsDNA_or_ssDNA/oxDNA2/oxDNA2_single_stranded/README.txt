
This is a simple example demonstrating how to build a short
single-stranded DNA polymer using the oxDNA2 force field with moltemplate.

This is based on the corresponding LAMMPS example located in this subdirectory:
examples/USER/cgdna/examples/oxDNA2/duplex1/
...which is distributed with LAMMPS.  (The conversion into MOLTEMPLATE format
was done automatically using the "ltemplify.py" file converter.)

Even though this is a single-stranded DNA molecule, in this example,
for its initial configuration, I decided to move and rotate each monomer (base)
into the same helical position it would be in if it were part of
ordinary (B-form) double-stranded DNA.  This is optional.
(See "polymer.lt" for details.)

WARNING: As of 2020-11-20, this example has not been tested or optimized.
         You should probably alter with the timestep, and langevin-settings
         to improve the simulation efficiency.

----
Note: If you get an error message (eg "Unknown pairstyle", "Unknown atom style")
      then you need to download the LAMMPS source code from git and
      compile LAMMPS with the following packages enabled:
             USER-CGDNA   ASPHERE
      See: https://lammps.sandia.gov/doc/Build_package.html
           https://lammps.sandia.gov/doc/Build_cmake.html
       (or https://lammps.sandia.gov/doc/Build_make.html)
----

Instructions on how to build LAMMPS input files and
run a short simulation are provided in other README files.

step 1)
README_setup.sh

step2)
README_run.sh
