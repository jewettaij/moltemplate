The purpose of this example is to test the density of water
constructed using the OPLSAA force-field.  (I think this is SPC water, not SPCE)

I just wanted some kind of sanity check to make sure we are converting
the OPLSAA parameters into moltemplate/LAMMPS format correctly.

The "TEST_density_estimate.txt" contains the results of that test.

-------- Instructions: ---------

More detailed instructions on how to build LAMMPS input files and
run a short simulation are provided in other README files.

step 1)
README_setup.sh

step 2)
README_run.sh


### Customizing atomic charges

In this example, atomic charge for OPLSAA atoms is determined by @atom type
(...according to a lookup table located at the beginning of the
"oplsaa2024.lt" file).
(Note: Any atomic charges listed in the "Data Atoms" section will be ignored.)
**These charges can be overridden.**
See the "README.md" file located in the parent directory
for instructions explaining how to customize atomic charge.
