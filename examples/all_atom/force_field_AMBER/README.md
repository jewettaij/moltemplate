This directory contains some examples of all-atom simulations using the
**AMBER** *(GAFF/GAFF2)* force field, prepared using moltemplate.

### Warning

This force field has not been extensively tested.  There is no gaurantee
that simulations prepared using moltemplate will exactly reproduce the
behavior of AmberTools/AMBER.

If you notice a problem with these examples, please report it.
(jewett.aij @ gmail.com)

### Limitations: *atomic charges*

Some force-fields (like COMPASS, and moltemplate's version of OPLSAA) include rules for assigning partial charges to atoms.  Most force fields, including AMBER GAFF and GAFF2 do not.  So GAFF and GAFF2 users will have to obtain atomic charges by some other means, probably by using 3rd-party tools.  (Alternatively, LAMMPS' [fix qeq/point](https://lammps.sandia.gov/doc/fix_qeq.html) feature can be used to assign partial charges, especially for simple molecules containing only C, H, O, N atoms.  If this fix is run infrequently, or run only once at the beginning of the simulation, then it should not slow the simulation down significantly.)

*(In these examples, I obtained partial charges from the OPLSAA
parameter file located [here](http://dasher.wustl.edu/tinker/distribution/params/oplsaa.prm).)*

### Improper angles

The style of improper interaction used by AMBER force fields depends on an
angle which depends on the order of the atoms surrounding the central atom.
When multiple atoms have the same type, this creates ambiguity in atom order.
Since there is no guarantee that moltemplate will choose the same atom order
as other molecule builders (such as AmberTools), this can lead to small
unavoidable discrepancies in energies and forces computed by LAMMPS and AMBER.
But their effect should be neglegible.

## Bloated lammps input scripts

By default, LAMMPS input scripts prepared using moltemplate contain the
entire contents of the GAFF or GAFF2 force-field, even when simulating small
systems with just a few atom types.

This is harmless, but if you want to get rid of this extra information,
follow the instructions in the "README_remove_irrelevant_info.sh" files.
