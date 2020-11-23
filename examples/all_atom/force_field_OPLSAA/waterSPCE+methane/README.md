Methane, water mixture
====================
This example contains a mixture of water(SPCE) and methane.  The methane molecules use OPLSAA force-field, but the water molecules do not.  The water molecules were initially arranged in a rectangular lattice.  The methane molecules were also arranged in a lattice, and were shifted to avoid overlap with the water molecules.  *(Alternatively, you can create a single lattice and specify the number of water and methane molecules you want in it using moltemplate's "new random([],[])" command, which is explained in the manual.  This gives you more control over the concentration of each ingredient.  You can also use PACKMOL to create random mixtures of molecules.)*


## Details

The methane molecules in this example use the OPLSAA force-field.  This means that the database of force-field parameters in the [oplsaa.lt file](../../../../moltemplate/force_fields/oplsaa.lt) will be used to generate angles, dihedrals, and impropers.  The [methane.lt file](moltemplate_files/methane.lt) contains these lines which refer to OPLSAA:
```
import "oplsaa.lt"
Methane inherits OPLSAA { ... }    # (see "methane.lt")
```
However the "SPCE" (water) molecules does NOT use a database to look up the force-field parameters for this tiny molecule.  Instead, the [spce.lt file](moltemplate_files/spce.lt) declares all of the angle interactions, atom properties and force-field parameters for water explicitly. (Consequently, it makes no mention of "OPLSAA".)

### Instructions

More detailed instructions on how to build LAMMPS input files and
run a short simulation are provided in other README files.

1) To build the files which LAMMPS needs, follow the instructions in:
[README_setup.sh](README_setup.sh)

2) To run LAMMPS with these files, follow these instructions:
[README_run.sh](README_run.sh)

(The instructions in "README_remove_irrelevant_info.sh" are optional.  *(If you notice a problem with this example, please [report it](../README.md).*)
