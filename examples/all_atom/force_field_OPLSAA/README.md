This directory contains some examples of all-atom simulations using the OPLSAA
force field.

### WARNING

This software is experimental and the force-fields and equilbration protocols
have not been tested carefully by me.  There is no gaurantee that simulations
prepared using moltemplate will reproduce the behavior of other MD codes.

If you notice a problem with these examples, please report it.
Peer-review is the only way to improve this software (or any software).
(jewett.aij @ gmail.com)

### Improper angles

The style of improper interaction used by OPLS force fields depends on an
angle which depends on the order of the atoms surrounding the central atom.
When multiple atoms have the same type, this creates ambiguity in atom order.
Since there is no guarantee that moltemplate will choose the same atom order
as other molecule builders (such as VMD), this can lead to small
unavoidable discrepancies in energies and forces computed by LAMMPS and NAMD.
But their effect should be neglegible.
*(Please let us know if this is not the case.)*

### Bloated lammps input scripts

By default, LAMMPS input scripts prepared using moltemplate contain the
entire contents of the OPLS force-field, even when simulating small
systems with just a few atom types.

This is harmless, but if you want to get rid of this extra information,
follow the instructions in the "README_remove_irrelevant_info.sh" files.
