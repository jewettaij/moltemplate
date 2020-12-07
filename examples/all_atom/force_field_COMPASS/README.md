This directory contains some examples of all-atom simulations using the COMPASS force field.

### WARNING

This software is experimental and the force-fields and equilbration protocols have not been tested carefully by me.  There is no gaurantee that simulations prepared using moltemplate will reproduce the behavior of other MD codes.

If you notice a problem with these examples, please report it. Peer-review is the only way to improve this software (or any software).

(jewett.aij @ gmail.com)

### Limitations

The moltemplate implementation of COMPASS currently relies on the same incomplete force-field file that "msi2lmp" uses ("compass_published.frc").  Unfortunately this means that many force field parameters and some atom types (such as sp2-carbons) have not (yet) been publicly released and are not available.

### Suggestions
Currently I recommend that users should run the "cleanup_moltemplate.sh" script after running "moltemplate.sh system.lt".  Then manually check that  the "system.in.settings" and "system.in.charges" files which remain make sense.  Specifically, you must check that the angle_coeff, dihedral_coeff, bond_coeff commands are not full of zeros (in places where they should not be zero.  This is another consequence of the fact that the .FRC files I mentioned above are incomplete.)  It's a good idea to also check that the charges in the "system.in.charges" file seem reasonable (ie. not all zeros).  (There is a list of warnings at the end of the "compass_published.lt" file.  You can check to see if any of the bonds in your system are covered by these warnings.)  Later on hopefully I'll add some automated way to warn users when these problems arise, but now you should check for them manually.


### Atomic charges

In most of the COMPASS examples,
the atomic charges are determined by their @atom types
and bond parters *(...according to the rules in the 
["compass_published.lt"](../../../moltemplate/force_fields/compass_published.lt) file)*.
*(Any atomic charges listed in the "Data Atoms" section of your molecules'
LT files will be ignored.)*
**These charges can be modified.**


### Customizing atomic charges in COMPASS molecules

#### Background information

LAMMPS provides two different methods to specify atomic charges:
1) **Specify charges in a DATA file** (eg "system.data").
*(This is the most popular way to specify atomic charges.
After running moltemplate.sh, the information in the "Data Atoms" section
is copied into the "Atoms" section of the "system.data" file created by
moltemplate.sh, and later read by LAMMPS.)*
2) **Specify charges using "set" commands.**
*(This is how the COMPASS atom charges are specified.
After running moltemplate, atom charge information in the
["compass_published.lt" file](../../../moltemplate/force_fields/compass_published.lt)
is copied into the "system.in.charges" file created by moltemplate.sh.
A LAMMPS input script file (eg. "run.in.nvt" or "run.in.npt")
is included with all of the COMPASS examples.  It tells LAMMPS to read
this "system.in.charges" file after reading the "system.data" file,
This overrides the atom charges from the "system.data" file.)*


#### How to customize atomic charge
*(without modifying "compass_published.lt")*

1) If you use use a 3rd-party program to calculate each atom's charge, you can
copy those charges into the "Data Atoms" section of your molecule's LT files.
To prevent LAMMPS from ignoring these charges, delete or comment-out the line
containing: **"include system.in.charges"** from your LAMMPS input script
(such as "run.in.min", "run.in.nvt", and "run.in.npt").
2) Alternatively, if you only want to override the charges of *some* of the
atoms in your molecules (and use default "compass_published.lt" charges for
the remaining atoms), then you can do this by adding an "In Charges" section
to your LT file and providing a list of custom charges for the \$atoms you
want to modify.  This is demonstrated in the
["graphene_nh2.lt"](../force_field_OPLSAA/functionalized_nanotubes_NH2/moltemplate_files/graphene_nh2.lt)
file located in
[this example](../force_field_OPLSAA/functionalized_nanotubes_NH2).
*(This example uses the OPLSAA force field, not the COMPASS force field,
but both force fields assign atomic charge according to atom type.
If you use this strategy, do not comment out "include system.in.charges"
from your "run.in\*" script files.)*
3) The discussion so far only applies to molecules that use the COMPASS
(or OPLSAA) force fields *(i.e. molecules whose definition begins with
"inherits COMPASS" or "inherits OPLSAA")*.
You can also mix molecules that use COMPASS or OPLSAA with other molecules
that don't.  (Those molecules will probably have their charges specified
in the "Data Atoms" section of their LT files.)
