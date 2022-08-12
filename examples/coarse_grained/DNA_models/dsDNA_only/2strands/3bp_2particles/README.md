3bp2p DNA model
===============

The moltemplate examples in this directory
all use the "3bp2p" model
(3 base-pairs per 2 beads)
which has an explicit double-helical shape.
This coarse-graind DNA model has been
[parameterized](simple_dna_example/moltemplate_files/deriving_force_field_parameters)
and
[tested](simple_dna_example#features).
It should be noted that there are much simpler
coarse-graind DNA models available with
similar mechanical properties.


### Prerequisites

LAMMPS must be compiled with the "MOLECULE" AND "EXTRA-MOLECULE"
packages enabled.  If LAMMPS generates the following error:
"dihedral_style spherical: Unknown dihedral style", then you must follow
[these instructions](https://lammps.sandia.gov/doc/Build_package.html),
and recompile LAMMPS.

