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

In most moltemplate examples, atomic charges (if present) are listed in
the 4th column of the "Data Atoms" section of each molecule's definition.
However the charges of atoms belonging to molecules which begin with
"inherits OPLSAA" is determined by their @atom types
*(according to a lookup table located at the beginning of the
["oplsaa.lt" file](../../../moltemplate/force_fields/oplsaa.lt) file)*.
**This can be overridden.**
See [here](../README.md#Customizing-atomic-charges)
for instructions how to customize atomic charges.
