This directory contains some examples of all-atom simulations using the OPLSAA force field.

### WARNING

There is no gaurantee that simulations prepared using moltemplate will reproduce the behavior of other MD codes.  If you notice a problem with these examples, please report it. Peer-review is the only way to improve this software (or any software).  (jewett.aij @ gmail.com)


### Atomic charges

In these examples, the charges in the "Data Atoms" section of each
molecule are ignored by LAMMPS.  That is because *(in this implementation of
OPLSAA)* the atomic charges are directly associated with their @atom types.
In other words, by default, two atoms with the same @atom type have the same
charges (regardless of the charges listed in the "Data Atoms" section).

**This can be overridden.**

If you use use a 3rd-party program to calculate each atomic charge, you can
copy those charges into the "Data Atoms" section of your molecule's LT files.
To prevent LAMMPS from ignoring the "Data Atoms" section, edit the
LAMMPS input script files (eg. "run.in.min", "run.in.nvt", "run.in.npt")
and delete or comment-out the line containing: **"include system.in.charges"**.
*(The "system.in.charges" file contains a series of commands will that assign atomic charges
according to their @atom types using the rules specified at the beginning of the
["oplsaa.lt" file](../../../moltemplate/force_fields/oplsaa.lt).)*

Alternatively, if you only want to override the charges of *some* of the atoms
in your molecule (and use default "oplsaa.lt" charges for the remaining atoms),
then you can do this by adding an "In Charges" section to your molecule
and provide a list of custom charges for the \$atoms you want to modify.
This is demonstrated in the ["graphene_nh2.lt"](functionalized_nanotubes_NH2/moltemplate_files/graphene_nh2.lt)
file located in [this example](functionalized_nanotubes_NH2).
*(In that example, the charge of one of the carbon atoms in the "Graphene_NH2"
 object was modified.
Note that if you do this, then do not comment out "include system.in.charges"
from all the "run.in\*" script files.)*



### Improper angles

The style of improper interaction used by OPLS force fields depends on an angle which depends on the order of the atoms surrounding the central atom. When multiple atoms have the same type, this creates ambiguity in atom order. Since there is no guarantee that moltemplate will choose the same atom order as other molecule builders (such as VMD), this can lead to small unavoidable discrepancies in energies and forces computed by LAMMPS and NAMD.  But their effect should be neglegible.
*(Please let us know if this is not the case.)*

### Bloated lammps input scripts

By default, LAMMPS input scripts prepared using moltemplate contain the entire contents of the OPLS force-field, even when simulating small systems with just a few atom types.

This is harmless, but if you want to get rid of this extra information, follow the instructions in the "README_remove_irrelevant_info.sh" files.
