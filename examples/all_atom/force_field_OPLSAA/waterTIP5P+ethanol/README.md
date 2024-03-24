TIP5P water, ethanol mixture
====================
This example of a mixture of rigid TIP5P water with an ordinary (flexible) molecule (ethanol) was contributed by Cheng Chen (Queen Mary University) on 2018-8-02.  The ethanol molecules use the OPLSAA force-field, but the water molecules do not.  The water molecules were initially arranged in a rectangular lattice.  The ethanol molecules were also arranged in a lattice, and were shifted to avoid overlap with the water molecules.  *(Alternatively, you can create a single lattice and specify the number of water and ethanol molecules you want in it using moltemplate's "new random([],[])" command, which is explained in [the manual](https://moltemplate.org/doc/moltemplate_manual.pdf#subsubsection.8.9.1).  This gives you more control over the concentration of each ingredient.  You can also use PACKMOL to create random mixtures of molecules.)*

The "fix rigid" command (or "fix rigid/nph/small") was used to integrate the equations of motion for the TIP5P molecule.  A combination of "fix nve" and "fix langevin" was used to integrate the equations of motion for the flexible ethanol molecule.

However there are multiple other ways to run these kinds of simulations.  A combination of "fix rigid/npt/small" and "fix nvt" might also work (without "fix langevin").  For a discussion of alternative ways to simulate mixtures of flexible and rigid molecules, see the [run.in.npt](run.in.npt) file.

### Instructions

More detailed instructions on how to build LAMMPS input files and
run a short simulation are provided in other README files.

1) To build the files which LAMMPS needs, follow the instructions in:
[README_setup.sh](README_setup.sh)

2) To run LAMMPS with these files, follow these instructions:
[README_run.sh](README_run.sh)

(The instructions in "README_remove_irrelevant_info.sh" are optional.)


### Customizing atomic charges

In this example, atomic charge for OPLSAA atoms is determined by @atom type
*(...according to a lookup table located at the beginning of the
["oplsaa.lt"](../../../moltemplate/force_fields/oplsaa.lt) file)*.
*(Any atomic charges listed in the "Data Atoms" section of your molecules'
LT files will be ignored.)*
**These charges can be overridden.**
See [here](../README.md#Customizing-atomic-charges-in-OPLSAA-molecules)
for instructions explaining how to customize atomic charge.
