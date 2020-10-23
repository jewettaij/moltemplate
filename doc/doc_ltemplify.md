ltemplify.py
===========


## Description

The "ltemplify.py" script is used to convert LAMMPS data files and
input scripts into a single MOLTEMPLATE ("**LT**") file.

*Typically*, the LT files generated
by "ltemplify.py" contain the definition of a single type of molecule
(or molecular complex) present in a data file (including geometry,
topology, force field information, and some of the the groups or fixes
that the molecule participates in.)  This way, moltemplate users later
on can build complicated simulations using this molecule as a building
block (perhaps along with other molecules).

Users can select the molecule (or molecules) they want
using the "-mol", "-id", or "-type" arguments,
and give the molecule a name using the "-name" argument.  (See below.)
The resulting LT file will include all of the information relevant
to that molecule including atom types, charges, coordinates,
bonded interactions, force field parameters, force field styles,
groups, and fixes that effect the molecule.
(Other information will be omitted.)

*However, by default*, "ltemplify.py" will copy all of the information
from these files into an LT file that
describes the entire system.
*(Later, when moltemplate.sh is run on that LT file,
  it will try to re-generate all of the original LAMMPS files.)*
Normally, this is not very useful.


Note: This program is both a stand-alone executable program (that can be run
from the terminal) and a python library.  The former is documented below.
*(The [python API is explained later](#Python-API).)*


### Typical Usage

```
ltemplify.py -name MoleculeName -mol MolID INPUT_SCRIPT DATA_FILE > FILE.lt
```
(...where MoleculeName is a string, MolID is an integer, INPUT_SCRIPT
and DATA_FILE are the names of a LAMMPS input script and a data file
containing the molecule of interest, and FILE.lt is the resulting
MOLTEMPLATE file created by ltemplify.py.)

***Note: Tiresome details to follow.***

### *First time readers should skip to the [examples section](#examples)*

--------------------------------------------------------

### Required arguments

"ltemplify.py" requires one argument: the name of a LAMMPS data file.
However (as shown in the example above), it also reads LAMMPS input scripts.
*(Note: If LAMMPS input scripts are included, then they must appear before
the DATA file in the argument list.  See examples below.)*
"ltemplify.py" also accepts many arguments to select the atoms belonging
to the molecule of interest and customize the output.


## Optional arguments

|Argument         | Explanation  |
|-----------------|--------------|
|-name *NAME*  | Specify the name of the molecule described in the output file. (The optional "inherits" keyword can be used to select a force field.) |
|-atomstyle "*style*" | Specify the name of the LAMMPS atom_style.  By default ltemplify.py will determine this from a comment following the "Atoms" line in the DATA file.  Use quotes to surround hybrid atom styles.  (eg. -atomstyle "hybrid bond ellipsoid").  If absent the "full" atom style is used by default.  |
|-columns "*column_list*" | As an alternative to -atomstyle, you can specify the name of each column in the Atoms section.  This list must be surrounded by quotes (eg -columns "id mol type charge x y z") |
|-mol "*molID_list*" | Select the molecule you want to appear in in the resulting LT file according to its molecule-ID number.  To select more than one molecule, supply a list of numbers surrounded by quotes. |
|-id "*id_list*" | Select atoms by their atom-IDs.  Use quotes to surround the list of atom ID numbers.  Unselected atoms will be omitted.  |
|-type "*type_list*" | Select atoms by their types. Use quotes to surround the list of atom type numbers. |
|-datacoeffs | Put force field information in the data file, not in an input script. (By default, force field parameters will be placed in the "In Settings"  section which will eventually be written to a LAMMPS input script.) |
|-ignore-comments  | Do not infer atom, bond, angle, dihedral, and improper type names from comments in the data file.  |
|-infer-comments  | Infer atom, bond, angle, dihedral, and improper type names from comments in the data file. |
|-ignore-coeffs  | Ignore all force field parameters (coeffs).  Omit from output file. (This is useful when using external force fields, such as OPLSAA.)|
|-ignore-angles  | Ignore angles, dihedrals and impropers.  Omit from output file. (This is useful when using external force fields, such as OPLSAA.)|
|-ignore-bond-types  | Ignore the 2nd column in the "Bonds" section of the LAMMPS data file, and create a "Bond List" section in the resulting MOLTEMPLATE LT file which omits the bond types.  (This is useful when using external force fields, such as OPLSAA.)| |
|-ignore-masses | Ignore all masses in the data file.  Omit from output file.  (This is useful when using external force fields, such as OPLSAA.)|
|-prepend-atom-type *STR* | prepend the string from the *STR* argument to the beginning of all atom type names. |
|-preamble *STR* | Print the string *STR* at the beginning of the file generated by ltemplify.py.  (Example: *'import "oplsaa.lt"'*)  This argument can be used multiple times if you want to prepend multiple strings at the beginning of the file. |

Examples showing typical argument usage are
[included below](#-Examples).

#### Default behavior

*Note that by default
(if the "-mol", "-id", or "-type" arguments are omitted),
"ltemplify.py" will copy all of the information
from the LAMMPS files into an LT file that describes the entire system.
Normally, this is not very useful.*


#### Details

All atoms, bonds, angles, dihedrals, and impropers and their associated
types will be converted to moltemplate "$" or "@" counter variables
(and the relevant portion of each file will be moved to sections
with the correct header names).
Coefficients, atom styles, and most force-field styles and
settings *should* also be included in the resulting .LT file.
ltemplify.py also understands simple group commands
(using "id", "molecule", or "type" styles)
and "fix shake" and "fix rigid"  (untested 2019-8-31).
However most other fixes, and complex group commands are not understood.
Those commands must be added to the resulting .LT file manually.
(More details [here](#-Known-bugs-and-limitations).)

### Fixes and Groups


ltemplify.py has *limited* support for "fix" and "group" commands,
including "fix shake", "fix rigid", and "fix poems".
Other fixes must be added manually to the file generated by ltemplify.py.
(Such as fix "restrain",
 "bond/create", "bond/break", "bond/react", "ttm", etc...)

ltemplify.py can understand simple (static) "group" commands, and will include them in the output file, if it can determine that they contain any relevant atoms.  (Fixes depending on irrelevant groups are also deleted.)


*Note: This feature has not been tested carefully.  So please review all of the group and fix commands generated by ltemplify.py to make sure they refer to the correct atoms.  And please report any bugs you find. (-Andrew 2018-8-28)*


### Automatic generation of atom, bond, angle, dihedral, improper names

By default ltemplify.py generates atom, bond, angle, dihedral, and improper,
type names and id names automatically.
This resultis in atoms with types like "@atom:type3", and IDs like
"$atom:type3_7" (I.e. the 7th atom of type 3.)

### Inferring atom type names from comments

*However,* ltemplify.py uses comments in the "Masses" section of
the LAMMPS DATA file (if present) to determine the name
of each atom type.
Consider the following excerpt from a hypothetical data file:
```
Masses

1 12.01  # c3
2 1.008  # h3
3 1.008  # ho
4 16.00  # oh
```
This means atoms of types 1, 2, 3, and 4 will be referred to as
"@atom:c3", "@atom:h3", "@atom:ho" and "@atom:oh",
respectively in the moltemplate (LT) file created by ltemplify.py.

#### Ignoring comments

The "*-ignore-comments*" argument will disable this behavior
and assign numeric names to the atom types in the usual way
(eg
"*@atom:type1*",
"*@atom:type2*",
"*@atom:type3*",
"*@atom:type4*").



### Bond, Angle, Dihedral, and Improper type names

Similarly, by default, bonds and angles are automatically
assigned to type names like "@bond:type4",
"@angle:type7".

*However, if comments appear* directly following the line in the
header file "*N* bond types", then these comments will be interpreted
as a list of bond type names (optionally preceded by an integer).
(The same is true of angle, dihedral, and improper type names.)
Consider this excerpt from a LAMMPS data file:
```
2 atom types
# c3
# h3

2 bond types
# CCethane
# c3_h3

2 angle types
# c3_c3_h3
# h3_c3_h3
```
In this example, bonds of type 1 and 2 will be referred to as
"@bond:CCethane" and "@bond:c3_h3"
in the moltemplate file, respectively.
Similarly, angles of type 1 and 2 will be referred to as
"@angle:c3_c3_h3" and "@angle:h3_c3_h3", respectively.
(As in the previos example, atoms of type 1 and 2 will be referred to as
"@atom:c3" and "@atom:h3" respectively.
You can specify atom type strings *either* here,
or in the Masses section.)


(As before, the "*-ignore-comments*" argument will disable this behavior.)

If you forget to add comments to the LAMMPS data file before running
*ltemplify.py*, you can always use a text-editor (or *sed*)
to manually find and replace all instances of "@atom:type1" with something
more meaningful, like "@atom:c3", for example.


## Force fields

Some data files contain a list of *angle, dihedral, or improper*
bonded interactions.  If so, then by default *ltemplify.py*
will include this information in the moltemplate (LT) file that it creates.
Sometimes, data files lack this information.

Either way, force fields
(including "OPLSAA", "GAFF2", and "COMPASS"),
contain rules for generating these interactions automatically.
Hence, users may intentionally wish to exclude this
information from the moltemplate files that ltemplify.py generates
when this information is contained in the force field they want to use.
(They can do this using the "-ignore-coeffs", "-ignore-angles",
 and "-ignore-bond-types" arguments explained below.)

### Using the inherits keyword to specify force fields

Moltemplate provides several different force fields to choose from
(such as OPLSAA, GAFF2, or COMPASS).  In addition, users can create their own custom force-fields.
To use these force fields, you must specify the one you want to use
using the *-name* argument with the *inherits* keyword
("**-name** "MOLECULE_NAME inherits FORCE_FIELD"")
For example:
```
ltemplify.py -name "Ethane inherits GAFF2" \
             -ignore-coeffs \
             ethane.data > ethane.lt
```
This will ask ltemplify.py to create a file defining
molecule named "Ethane".
Later when moltemplate is used to read this file, the "GAFF2"
force field will be used to generate angles, dihedrals and impropers,
and lookup their force field parameters.

In addition, after ltemplify.py is finished, the user must manually insert the following line *at the beginning* of the file that ltemplify.py created.  For example:
```
import "gaff2.lt"  #<-- define the GAFF2 force field so we can use it later

# --- the text below was generated by ltemplify.py ---
Ethane inherits GAFF2 {
  ...
}
```
ltemplify.py does not do this for you.  A list of available force fields can be found in the "moltemplate/force_fields/" directory distributed with moltemplate on github.



### -ignore-coeffs

The optional "*-ignore-coeffs*" argument will
force ltemplify.py to ignore the  force field parameters
that it encountered in the user's input script or DATA file.
The resulting LT file will omit this information.
If you plan to use a force field with this molecule, then this information
will be present in the force field you are using, so there's no need
to include it in the resulting LT file you are creating now.
(Later when you run moltemplate.sh on the LT file that ltemplify.py
created, it will use the force field to lookup these force field parameters.)

### -ignore-angles

If the original DATA file has "Angles", "Dihedrals", or "Impropers",
you can use the "*-ignore-angles*" argument if you want to force
ltemplify.py to ignore/remove those interactions from the LT
file which ltemplify creates.
(Doing that will allow the force field rules to take precedence
later when we run moltemplate.sh on that file.)

### -ignore-bond-types

Similarly, when using force-fields, you only need to specify a
list of *which pairs of atoms* are bonded together.
The force-field will determine the type and
properties of each bond (eg, equilibrium rest length, stiffness, etc...)
according to atom type names and the force field rules.

To do that, you must force *ltemplify.py* to ignore the existing
bond type information present in your data file using the
"*-ignore-bond-types*" argument.
This will force ltemplify.py to ignore the bond types in the
(2nd column of the) "Bonds" section of the LAMMPS data file that you provided.
In this way, the bond type can be determined later by moltemplate.sh
in a way which is consistent with the force field you selected.


See [below](#-Examples-using-force-fields) for examples.



#### Warning

*ltemplify.py is experimental software.*
*Unlike moltemplate.sh*, the *ltemplify.py* script has
limited understanding of all of the features available in LAMMPS.
Please review the resulting ".LT" file and check for errors.
(If necessary, convert any remaining
atom, bond, angle, dihedral, or improper id or type numbers to the
corresponding \$ or @ variables.)
Some exotic pair styles which have their own special syntax
are not understood.
These interactions must be converted to moltemplate format manually.
Support for "group" and "fix" commands is also
[experimental](#-Fixes-and-Groups)
Please report errors in the behavior of ltemplify.py.



## Examples


#### Example 1


```
ltemplify.py -name Ethane -molid "1" FILE.in FILE.data > ethane.lt
```

This example creates a new file ("ethane.lt")
containing a new type of molecule (named "Ethane"),
consisting of all the atoms whose molecule-ID number equals 1.
*(Presumabely, the first molecule in FILE.data is an ethane molecule.)*

ltemplify.py reads the atom coordinates and bonded interactions from FILE.data.
Other information relevant to that molecule (including the atom_style,
force-field styles and parameters, groups and fixes)
are read from "FILE.in" (which is presumabely a LAMMPS input script file).

*(NOTE: Again, it is not necessary to include a LAMMPS input script in
        the argument list.  However important information is typically
        contained in LAMMPS input script files, so if you have one, including
        it is recommended.  However a data file is enough.)*

Note: Selecting atoms by molecule-ID only works if you are using
one of the "molecular" atom\_styles (such as "atom\_style full").
If you are using a different atom\_style
(such as "atom\_style angle" or "atom\_style bond"),
you can select the atoms you want either by type or by id number.
(See below.)



#### Example 2

Sometimes, the information describing your molecule will divided
into multiple lammps input scripts.
(For example, one input script may contain various *style* commands.
 The next input script may contain *coeff* commands.)
In that case, these input scripts should appear
in the argument list *before the data file*,
and in the order in which they are read by LAMMPS.
```
ltemplify.py -name Ethane -molid "1" FILE1.in FILE2.in FILE.data > ethane.lt
```



#### Example 3

```
ltemplify.py -name Ethane -id "13 14 15 61*69" FILE.in FILE.data > ethane.lt
```

    In this example, only atoms whose ids are
    13, 14, 15, and 61 through 69 are included.


#### Example 4

```
ltemplify.py -name Ethane -type "1 2 3" FILE.in FILE.data > ethane.lt
```

    In this example, only atoms whose type is 1, 2, or 3 are included.




#### Example 5

```
ltemplify.py -name EntireSystem FILE.in FILE.data > entire_system.lt
```

This creates a template for a new molecule object (named "EntireSystem"),
consisting of ***all*** the atoms in the lammps files you included,
and saves this data in a single LT file ("entire\_system.lt").
This file can be used with moltemplate.sh (and/or ttree.py) to
define large systems containing this molecule.

Note: Again, the input scripts ("FILE.in" in this example) should appear 
      before the data file ("FILE.data") in the argument list.




### Examples using force fields:

You can also use *ltemplify.py* to create molecules that use 3rd-party
force fields such as OPLSAA, GAFF2, COMPASS, ....

#### Example 6

This example demonstrates how to build a molecule
using the "GAFF2" force field.
The following example extracts molecule 1 from "FILE.in" and "FILE.data".

```
# This example creates a new file, "ethane.lt", which will contain the
# instructions for building a "Ethane" molecule using "GAFF2". First
# specify which file contains the definition of the "GAFF2" force field:

echo "import gaff2.lt"  >  ethane.lt

# Then use ltemplify.py to extract information from FILE.in, FILE.data

ltemplify.py -name "Ethane inherits GAFF2" \
             -molid "1" \
             -ignore-angles -ignore-bond-types -ignore-coeffs \
             FILE.in FILE.data >> ethane.lt

# Note: if you want to build a simulation containing these molecules,
# you will have to create a "system.lt" file which refers to "ethane.lt"
# and then run moltemplate.sh on this file.
```

As mentioned earlier, comments in "file.data" will determine the name
of each atom type and *should match atom type names in the force field*.

In this example, the angle, dihedral, improper, and bond-type
information is stripped from the original file.data
(and will be generated later according the the rules
defined in the "GAFF'2' force field).
The name of the molecule ("Ethane inherits GAFF2") includes a reference
to the force field ("GAFF2") which will be used to lookup this information.
(Note: The "GAFF2" force field parameters are typically defined in a file
       named "gaff2.lt".  Hence in this example we used "echo"
       to insert a link to "gaff2.lt" at the beginning of the
       "ethane.lt" file so that moltemplate.sh will know where to find them.
       Alternatively, this could be done manually by the user.)


## Known bugs and limitations

#### Exotic styles are not supported

ltemplify.py does ***not*** understand the syntax of
exotic many-body pair\_styles such as tersoff, sw, meam, reax, dpd, edip,
dipole, lubricate, hbond/dreiding
(even though these styles are supported by moltemplate).
After running ltemplify.py, the user must manually edit the resulting ".lt"
files.  For example: ltemplify.py will not understand wildcard characters
("*" characters)
which typically appear in the "pair\_coeff" commands or "Pair Coeffs"
section when using these many-body pair styles.
You will have to remove the extra lines automatically generated by ltemplify.py
and put the wildcard characters back (eg "pair\_coeff * * ...") manually.
(Later the user may need to run moltemplate using the appropriate "-a"
 command line args to make sure the various atom types are assigned
 to the correct numbers. This is usually needed in order to keep them
 consistent with the order of parameters in the corresponding pair style's
 input files.  In moltemplate you can manually assign atom types to numbers
 using the *-a* argument.)
In addition, auxiliary atom types (such as the "hydrogen" atom type
required by hbond/dreiding) will not even be parsed.
If you are using the "hbond/dreiding" pair style, you will
have to manually specify the atom type for the hydrogen-atom mediator
in every "pair\_coeff" command after running ltemplify.py


#### Wildcard characters ("*") expansion

Moltemplate is often confused whenever wildcard characters ("*" characters)
appear inside any of the the "coeff" commands
(or "Coeff" sections of the data file).
So ltemplify.py attempts to remove these characters and expand these commands,
generating multiple lines of output, and listing each atom type explicitly.
(This is also done for bond types, angle types, dihedral types,
 and improper types.)
This may not be what you want.
(For example, this can be a problem if you are using a many-body pair style
which requires you to specify "* *" for the atom types, such as
*tersoff*, *eam*, or *sw*.)


## Python API

It is possible to access the functionality of *ltemplify.py* from 
within python  (*without* invoking ltemplify.py using subprocess.run(), 
os.system(), or writing files to the file system).
To do that, you can create python strings containing the 
contents of the LAMMPS data file (and input scripts) and use 
Ltemplify.Convert() to convert them into MOLTEMPLATE (.LT) format.  
 
However, (as you can probably tell) making this possible within python 
was an afterthought.  Currently, the recommended way to do this is to pass the 
same command line arguments (described above) to the constructor of the
*Ltemplify* object.  *Then* invoke the *Convert()* function.

*(Alternatively, you can also edit the data members of the Ltemplify object
 directly after it is created, instead of using the command-line arguments.
 But doing it that way exposes you to the horrifically messy contents
 of the Ltemplify object.  Perhaps in the future, I will clean it up.)*

```python
class Ltemplify(object):
    def __init__(self, argv):
        """
        The constructor requires a list of command line arguments to 
        figure out the format of the output file we will generate.
        This meaning of these arguments is explained above.
        Note: You can either specify the input scripts and data files in the
        argument list, OR specify them later by passing them as arguments
        to the Convert() member function.  (The second approach is preferred.)
        """
    def Convert(self,
                out_file,
                input_data_file=None,
                input_script_files=None):
        """
        Converts a data file (and, optionally, one or more input scripts)
        into a new file ("out_file") which is in MOLTEMPLATE (.LT) format.
        The arguments can be either filenames or StringIO objects.
        The "input_script_file" argument can be a single string or StringIO
        object, or a list of such objects.
        """
```

## Usage example inside python

The goal of this example is to demonstrate how to invoke the features of
*ltemplify.py* from within python using strings instead of files.

So I created long strings containing the contents of a LAMMPS data file,
and several input scripts.  Then I convert these files into StringIO objects.
Then I use Ltemplify.Convert() to convert them into MOLTEMPLATE (.LT) format.

*(My apologies for this example being so long.)*


```python
data_file_contents = \
"""
LAMMPS Description

     6  atoms
     5  bonds

     2  atom types
     # c2
     # hc

     2  bond types
     # C-C
     # C-H
     
  0.0 48.00 xlo xhi
  0.0 48.00 ylo yhi
  0.0 48.00 zlo zhi

Masses

1 12.011  # c2
2 1.008   # hc

Atoms  # full

1 1 1 0.00 -0.6695 0.0 0.0
2 1 1 0.00 0.6695 0.0 0.0
3 1 2 0.00 -1.23422 -0.85446 0.0
4 1 2 0.00 -1.23422 0.85446 0.0
5 1 2 0.00 1.23422 -0.85446 0.0
6 1 2 0.00 1.23422 0.85446 0.0

Bonds

1 2 1 2
2 1 1 3
3 1 1 4
4 1 2 5
5 1 2 6

"""

input_script_file1_contents = \
"""
  atom_style full
  units real
  bond_style harmonic
  pair_style lj/cut/coul/long 10.0 10.0
  pair_modify mix geometric
  special_bonds lj/coul 0.0 0.0 0.5
  kspace_style pppm 0.0001
"""

input_script_file2_contents = \
"""
  pair_coeff 1 1 0.076 3.55
  pair_coeff 2 2 0.03 2.42
  bond_coeff 1 340.0 1.08
  bond_coeff 2 549.0 1.34
  set type 1 charge -0.23
  set type 2 charge 0.115
  timestep   2.0
  dump   1 all custom 5000 traj_nvt.lammpstrj id mol type x y z ix iy iz
  fix    fxnvt all nvt temp 300.0 300.0 500.0 tchain 1
  run    200000
"""

import io
data_file          = io.StringIO(data_file_contents)
input_script_file1 = io.StringIO(input_script_file1_contents)
input_script_file2 = io.StringIO(input_script_file2_contents)

input_script_files = [input_script_file1, input_script_file2]

args=['-atomstyle', 'full',
      '-name','Ethylene inherits GAFF2',
      '-mol', '1',
      '-ignore-angles',
      '-ignore-bond-types',
      '-ignore-coeffs']

# Create an Ltemplify object with these settings:

import moltemplate
ltmp = moltemplate.ltemplify.Ltemplify(args)

output_file = io.StringIO()
# Now convert this to a file in MOLTEMPLATE format
ltmp.Convert(output_file, data_file, input_script_files)

output_file.seek(0)
output_file_contents = output_file.read()

# ("output_file_contents" is a string containing the contents of the
#  file in MOLTEMPLATE (.LT) format.)

```

