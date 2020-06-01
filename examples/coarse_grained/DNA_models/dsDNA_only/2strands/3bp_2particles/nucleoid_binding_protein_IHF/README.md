IHF nucleoid binding protein (NAP) example
====================

This example explores how a coarse-grained protein that causes 180-degree
bends in DNA might effect the conformation of (bacterial) chromosomes.
The coarse grained protein was inspired by the "IHF" protein which promiscuously
binds to DNA, at numbers approaching 30000 per cell in rapidly dividing E.coli.

The simulation begins with several hundred IHF proteins attached to a small
interval of DNA (177kbp) which has been stretched to it's full length.
A tension in the DNA is maintained (by running the simulation at negative
pressure in the X direction).  Gradually the tension in the DNA is reduced
to 0.  During this time, twist motors (distributed evenly along the DNA)
apply torsional torque to the polymer to encourage it to form supercoils.

We can vary number of IHF protiens and the torque in the twist motors,
and the pressure in the X direction to see how the effect of IHF protein
on the chromosome conformation varies with concentation, supercoil density,
and pressure.

Note: This simulation protocol assumes that IHF binds to DNA irreversibly.
More general interactions between NAPs and DNA may be easier to simulate
(at this low resolution) with the aid of LAMMPS features like
fix bond/react or other fixes which have similar features.

##    Prerequisites

LAMMPS must be compiled with the "MOLECULE" AND "USER-MISC" packages enabled.
(https://lammps.sandia.gov/doc/Build_package.html)

It also requires that the "fix twist" feature has been enabled in LAMMPS.
(As of 2019-5-05, you must download "fix_twist.cpp" and "fix_twist.h" from
 https://github.com/jewettaij/lammps/tree/fix_twist/src/USER-MISC
 and copy those 2 files into the "src/" subdirectory of our LAMMPS folder,
 and re-compile LAMMPS.  Hopefully in the future this won't be necessary.)

After enabling the packages you need (and, if necessary copying the
"fix_twist.cpp" and "fix_twist.h" files), you must (re)compile LAMMPS
to enable the features that this example uses.

If, when running LAMMPS, you receive this error message
"Unknown dihedral style", "Unknown fix", or something similar,
it means you did not successfully follow the instructions above.

##    WARNING

These files (originally uploaded on 2019-12-10) contain many comments
(beginning with "#") which are probably misleading and no longer relevant.
I will try to clean up these files over time.
Please let me know if anything doesn't work.

##    Instructions

Instructions on how to build LAMMPS input files and 
run a short simulation are provided below

The following file contain instructions explaining how to generate
the curve that you want the DNA polymer to follow.
(You can also run it as an executable.)

   ./STEP_1_generate_initial_path.sh

The next file explains how to convert this curve into a moltemplate file, and
how to run moltemplate on that file. (You can also run it as an executable.)

   ./STEP_2_generate_LAMMPS_files.sh

Finally, to run the LAMMPS simulation follow the instructions in this file:
STEP_3_run_sim.sh
You will have to edit the file to specify the name of the LAMMPS binary
you are using (for example, "lmp_ubuntu"), and the number of processors.

