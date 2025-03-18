This directory contains some examples of all-atom simulations using the OPLSAA force field.

### WARNING

There is no guarantee that simulations prepared using moltemplate will reproduce the behavior of other MD codes.  If you notice a problem with these examples, please report it. Peer-review is the only way to improve this software (or any software).  (jewett.aij @ gmail.com)



-------------------------

## Atomic charges

In most of the OPLSAA examples,
the atomic charges are determined by their @atom types
*(...according to a lookup table located at the beginning of the
["oplsaa2024.lt"](../../../moltemplate/force_fields/oplsaa2024.lt) file)*.
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
["oplsaa2024.lt" file](../../../moltemplate/force_fields/oplsaa2024.lt)
is copied into the "system.in.charges" file created by moltemplate.sh.
A LAMMPS input script file (eg. "run.in.nvt" or "run.in.npt")
is included with all of the OPLSAA examples.  It tells LAMMPS to read
this "system.in.charges" file after reading the "system.data" file,
This overrides the atom charges from the "system.data" file.)*


#### How to customize atomic charge
*(without modifying "oplsaa2024.lt")*

1) If you use use a 3rd-party program to calculate each atom's charge, you can
copy those charges into the "Data Atoms" section of your molecule's LT files.
To prevent LAMMPS from ignoring these charges, delete or comment-out the line
containing: **"include system.in.charges"** from your LAMMPS input script
(such as "run.in.min", "run.in.nvt", and "run.in.npt").
2) Alternatively, if you only want to override the charges of *some* of the
atoms in your molecules (and use default "oplsaa2024.lt" charges for the remaining
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



-------------------------

## Minor issue: Bloated lammps input scripts

By default, LAMMPS input scripts prepared using moltemplate contain
the entire contents of the OPLSAA force-field, even when simulating
small systems with just a few atom types.
This is harmless, but if you want to get rid of this extra information,
follow the instructions in the "README_remove_irrelevant_info.sh" files.


-------------------------

## *OPTIONAL*
## Customizing dihedrals, angles, and bonds

The OPLSAA force field contains many alternative parameter choices for
dihedral, angle, and bond interactions.
It's not always possible to for determine the optimal choice of dihedral angle
parameters from the @atom types alone.  By default, moltemplate hides this
issue, and will *attempt to make a reasonable guess*, chosing the most generic
version of the interaction between those atom types.
*(The same is true with most other molecule builder programs.)*

***Most of the time, this is fine.***
At the beginning when you are trying to get your simulation to run,
don't worry about choosing the optimal dihedral, angle, or bond parameters.
The default choice is often good enough.

*Eventually, if you want to improve the accuracy of your simulation somewhat,
then you can detect these ambiguous duplicate dihedrals and angles.
Then you can override the default choice by following the (somewhat laborious)
procedure below.*

Again, moltemplate hides this issue by default.
To be informed when moltemplate detects multiple ambiguous dihedrals
(or angles), you must run moltemplate.sh with the optional
`-report-duplicates bytype __` arguments.
For example:
```
moltemplate.sh  system.lt  -report-duplicates bytype __
```
- If you see a file named "warning_duplicate_dihedrals.txt" after running
moltemplate.sh, then moltemplate found multiple plausible
dihedral interactions between the same set of atoms in your molecules.
If you see this file, then read the first few warning messages in that file.
*(This file is typically long and redundant, since there are typically many
copies of the same molecules in a simulation.)*
Then modify your .lt files accordingly.
*(This requires adding a custom "Data Dihedrals" section where you specify
the version of the @dihedral that you want to use for those atoms.
Several example .lt files demonstrate how to do that in detail, including
[butane.lt](./butane/moltemplate_files/butane.lt) and
[benzoic_acid_optimizations.lt](./benzene+benzoic_acid/moltemplate_files/optimized_version_using_custom_dihedrals/benzoic_acid_optimizations.lt).)*
- Each time you add a line to your "Data Dihedrals" section, the corresponding
warning(s) in the "warning_duplicate_dihedrals.txt" file will be removed.
So read the first warning, add a corresponding line to the "Data Dihedrals"
section to address that warning, and run moltemplate.sh again.  Repeat this
until the "warning_duplicate_dihedrals.txt" file is no longer being generated.
- If you see a file named "warning_duplicate_angles.txt"
or "warning_duplicate_bonds.txt", then follow the same procedure.
Read those warning messages and add a custom "Data Angles" or "Data Bonds"
section to your .lt files to override the default choice of
@angle or @bond type.

