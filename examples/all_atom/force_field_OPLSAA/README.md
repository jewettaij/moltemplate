This directory contains some examples of all-atom simulations using the OPLSAA force field.

### WARNING

There is no gaurantee that simulations prepared using moltemplate will reproduce the behavior of other MD codes.  If you notice a problem with these examples, please report it. Peer-review is the only way to improve this software (or any software).  (jewett.aij @ gmail.com)


### Atomic charges

In most of the OPLSAA examples,
the atomic charges are determined by their @atom types
*(...according to a lookup table located at the beginning of the
["oplsaa.lt"](../../../moltemplate/force_fields/oplsaa.lt) file)*.
*(Any atomic charges listed in the "Data Atoms" section of your molecule's
LT files will be ignored.)*
**These charges can be overridden.**


### Customizing atomic charges in OPLSAA molecules

#### Background information

LAMMPS provides two different methods to specify atomic charges:
1) **Specify charges in a DATA file** (eg "system.data").
*(This is the most popular way to specify atomic charges.
After running moltemplate.sh, the information in the "Data Atoms" section
is copied into the "Atoms" section of the "system.data" file created by
moltemplate.sh, and later read by LAMMPS.)*
2) **Specify charges using "set" commands.**
*(This is how the OPLSAA atom charges are specified.
After running moltemplate, atom charge information in the
["oplsaa.lt" file](../../../moltemplate/force_fields/oplsaa.lt)
is copied into the "system.in.charges" file created by moltemplate.sh.
A LAMMPS input script file (eg. "run.in.nvt" or "run.in.npt")
is included with all of the OPLSAA examples.  It tells LAMMPS to read
this "system.in.charges" file after reading the "system.data" file,
This overrides the atom charges from the "system.data" file.)*


#### How to customize atomic charge
*(without modifying "oplsaa.lt")*

1) If you use use a 3rd-party program to calculate each atom's charge, you can
copy those charges into the "Data Atoms" section of your molecule's LT files.
To prevent LAMMPS from ignoring these charges, delete or comment-out the line
containing: **"include system.in.charges"** from your LAMMPS input script
(such as "run.in.min", "run.in.nvt", and "run.in.npt").
2) Alternatively, if you only want to override the charges of *some* of the
atoms in your molecules (and use default "oplsaa.lt" charges for the remaining
atoms), then you can do this by adding an "In Charges" section to your LT file
and providing a list of custom charges for the \$atoms you want to modify.
This is demonstrated in the ["graphene_nh2.lt"](functionalized_nanotubes_NH2/moltemplate_files/graphene_nh2.lt)
file located in [this example](functionalized_nanotubes_NH2).
*(In that example, the charge of one of the carbon atoms in the "Graphene_NH2"
object was modified.  If you use this strategy, do not comment out
"include system.in.charges" from your "run.in\*" script files.)*
3) The discussion so far only applies to molecules that use the OPLSAA force
field *(i.e. molecules whose definition begins with "inherits OPLSAA")*.
You can also mix molecules that use OPLSAA with other molecules
that don't.  In the [waterSPCE+methane](waterSPCE+methane) example,
the SPC/E water molecules do not use OPLSAA.
Hence, their atomic charges are located in the "Data Atoms" section
of the [spce.lt](waterSPCE+methane/moltemplate_files/spce.lt) file.
*(The same is true of most of the carbon atoms in the
[carbon nanotube](functionalized_nanotubes_NH2) example.)*



### Improper angles

The style of improper interaction used by OPLS force fields depends on an angle which depends on the order of the atoms surrounding the central atom. When multiple atoms have the same type, this creates ambiguity in atom order. Since there is no guarantee that moltemplate will choose the same atom order as other molecule builders (such as VMD), this can lead to small unavoidable discrepancies in energies and forces computed by LAMMPS and NAMD.  But their effect should be neglegible.
*(Please let us know if this is not the case.)*

### Bloated lammps input scripts

By default, LAMMPS input scripts prepared using moltemplate contain the entire contents of the OPLS force-field, even when simulating small systems with just a few atom types.

This is harmless, but if you want to get rid of this extra information, follow the instructions in the "README_remove_irrelevant_info.sh" files.
