This directory contains some examples of all-atom simulations using the OPLSUA
force field.

*Note:* Confusingly, the parameters and atom type definitions for the OPLSUA force-field are contained in the [oplsaa.lt](../../../moltemplate/force_fields/oplsaa.lt) file.  The first 56 atom types defined in that file (specifically *@atom:1*...*@atom:56*) correspond to OPLSUA united-atom particles.  To use the OPLSUA force field, select atoms for your molecule entirely from that list.

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
