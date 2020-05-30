genpoly_modify_lt.py
===========


## Description

Modify an existing polymer by adding (or overriding) bonds, or bonded
interactions, changing atom types, or adding constraint fixes
at *selected monomers* within an existing polymer created with
[genpoly_lt.py](doc_genpoly_lt.md).
Like genpoly_lt.py, this program creates (or augments) a moltemplate LT file
containing the definition of a polymer.  You will have to use moltemplate.sh
to convert this polymer(s) into a system that can be simulated in LAMMPS.
(This program does not communicate with LAMMPS directly.)

Unlike *genpoly_lt.py*, these modifications are typically made to
only a subset of the monomers in the polymer.
This program is a convenient way to adding the same modification
to a polymer in many places along its length.
The user typically supplies a file specifying the locations along
the polymer where you would like the modification to occur.
(You can also distribute them automatically at periodic or randomly
 spaced locations.)

When you run this program, the same modification is made
(typically multipe times) along the length of your polymer.
However if you need to make multiple different *kinds* of modifications
to the polymer, you can invoke this program multiple times at different
locations, concatenating the results each time to the end of the LT file
that describes your polymer


### Purpose

This program was originally written to help build complicated chromosomal
DNA examples with many modifications and bonded proteins distributed along
the length of the DNA.  I don't know if it is useful enough to be of general
interest outside of the community of people who simulate coarse grained DNA.
However this program can be used to modify any polymer (not just DNA).


## Usage:

```
      [-polymer-name pname] \\
      [-length num_monomers] \\
      [-locations filename] \\
      [-locations-periodic num_mods offset] \\
      [-locations-random num_mods seed] \\
      [-width mod_width] \\
      [-bond btype a1 a2 i1 i2] \\
      [-angle    atype a1 a2 a3 i1 i2 i3] \\
      [-dihedral dtype a1 a2 a3 a4 i1 i2 i3 i4] \\
      [-improper itype a1 a2 a3 a4 i1 i2 i3 i4] \\
      [-set-atoms M filename attribute a1 ... aM i1 ... iM A1 ... Am] \\
      [-fix-nbody N filename fixname fixID group keyword a1 ... aN i1 ... iN params] \\
      [-circular yes/no/connected] \\
      >> polymer.lt
```

## Arguments

```
    -polymer-name name
             Name of the moltemplate object that you created earlier.
             (If it was created using genpoly_lt.py, then it should match
              the argument supplied to its "-name" argument.)
                   
    -length num_monomers
             Specify the total number of monomers in the polymer.  (This is
             not optional.  This program has no way of knowing the length of
             the polymer you created earlier.)

    -locations locations.txt
             Supply a file containing a list of integers, starting at 0,
             indicating the locations along the polymer where you would like
             the modifications to the polymer to be made.

    -locations-periodic num_mods offset
             Alternatively, you can distribute these modifications at num_mods
             evenly spaced intervals along the length of the polymer, starting
             at position i=offset.  (Setting offset=0 means that the first
             monomer will be modified.)  The length of the polymer need not
             be a multiple of num_mods.

    -locations-random num_mods seed
             Alternatively, you can distribute these modifications randomly
             along the length of the polymer.  In addition to the number of
             modifications, you must also specify a seed for the random number
             generator. (Any integer will do.)  Care is taken to avoid overlap
             with each other.  You can control the minimum spacing between
             modifications using the "-width" argument.

    -width mod_width
             Specify the minimum "width" of a modification generated using
             the "-locations-random" argument.  Equivalently, this is (1+) the
             minimum space allowed between randomly generated modifications.
             (If left unspecified, it is 1 by default.  This means that the
             modification can be placed at successive monomers along the chain.)
             If you plan to add multiple different kinds of modifications,
             you can use the -widths argument instead.  (See examples below.)

    -bond btype a1 a2 i1 i2
             Add a bond interaction between atoms a1 and a2 in monomers, located
             at position i1 and i2 relative to the location where the
             modification occured. (The a1 and a2 arguments are strings
             indicating the atom-ID within the monomer
             omit the @bond: and $atom: prefixes in moltemplate variables)
             The i1, i2 arguments are integer indices indicating the (relative)
             location of monomer that each atom belongs to
                0 corresponds to the monomer specified in the "-locations" file
                1 corresponds to the next monomer
                2 corresponds to the following monomer, etc...
	     (For circular polymers, the indices will be wrapped appropriately.)
             These indices are added to the location where you want the
             modification to be made (as specified, for example, in the file
             supplied to the "-locations" argument, or generated automatically).
             Multiple bonds can be added by having "-bond btype a1 a2 i1 i2"
             appear several times in the argument list.
             Note that by default, when multiple bonds exist between the same
             pair of atoms the later bond overrides the earlier bond (unless
             -overlay-bonds is used).  Consequently, this provides a mechanism
             to disable or weaken existing bonds without actually deleting them.
             (You can delete them too using moltemplate's "delete" command,
              but this program does not generate the text for those commands.)

    -angle atype a1 a2 a3 i1 i2 i3
             Add a 3-body angle interaction between atoms a1 a2 a3 in monomers
             i1 i2 and i3.  (The atype a1, a2, a3 arguments are strings
             containing moltemplate variable names. The standard moltemplate
             prefixes "$angle:", "@angle:", and "$atom:" should be omitted.
             The i1, i2, i3 arguments are integer indices indicating the monomer
             that each atom belongs to.
                0 corresponds to the current monomer
                1 corresponds to the next monomer
                2 corresponds to the following monomer, etc...
	     (For circular polymers, the indices will be wrapped appropriately.)
             Multiple angles per monomer can be added by having:
             "-angle aname atype a1 a2 a3 i1 i2 i3"
             appear several times in the argument list with different parameters
             Note that by default, when multiple angle interactions exist
             between the same triplet of atoms, the later angle overrides the
             earlier angle (unless -overlay-angles is used).  Consequently,
             this provides a mechanism to disable or weaken existing angles
             without actually deleting them.  (You can delete them manually
             using moltemplate's "delete" command, but not using this program.)

    -dihedral dtype a1 a2 a3 a4 i1 i2 i3 i4
             Add a 4-body dihedral interaction between atoms a1 a2 a3 a4 in
             monomers i1 i2 and i3.  (The dname dtype a1, a2, a3, a4, arguments
             are strings containing moltemplate variable names. The moltemplate
             prefixes "$dihedral:", "@dihedral:", and "$atom:" should be omitted
             The i1, i2, i3, i4 arguments are integer indices indicating the
             monomer that each atom belongs to.  (See explanation above.)
             Multiple dihedrals per monomer can be added by having:
             "-dihedral dname dtype a1 a2 a3 a4 i1 i2 i3 i4"
             appear several times in the argument list with different parameters
             Note that by default, when multiple dihedral interactions exist
             between the same quartet of atoms, the later dihedral overrides the
             earlier dihedral (unless -overlay-dihedrals is used). Consequently,
             this provides a mechanism to disable or weaken existing dihedrals
             without actually deleting them.  (You can delete them manually
             using moltemplate's "delete" command, but not using this program.)

    -improper itype a1 a2 a3 a4 i1 i2 i3 i4
             Add a 4-body improper interaction between atoms a1 a2 a3 a4 in
             monomers i1 i2 and i3.  (The iname itype a1, a2, a3, a4, arguments
             are strings containing moltemplate variable names. The moltemplate
             prefixes "$improper:", "@improper:", and "$atom:" should be omitted
             The i1, i2, i3, i4 arguments are integer indices indicating the
             that each atom belongs to.  (See explanation above.)
             Multiple impropers per monomer can be added by having:
             "-improper iname itype a1 a2 a3 a4 i1 i2 i3 i4"
             appear several times in the argument list with different parameters
             Note that by default, when multiple improper interactions exist
             between the same quartet of atoms, the later improper overrides the
             earlier improper (unless -overlay-impropers is used). Consequently,
             this provides a mechanism to disable or weaken existing impropers
             without actually deleting them.  (You can delete them manually
             using moltemplate's "delete" command, but not using this program.)

    -set-atoms M filename attribute a1 ... aM i1 ... iM A1 ... Am
          Generate a file containing LAMMPS commands which will modify
          attributes of M atoms at each location specified by the user.
          This is done using the LAMMPS' "set" command.
          https://lammps.sandia.gov/doc/set.html
          "M" is the number of atoms you want to modify at each location.
          "filename" is the name of a file where the LAMMPS commands will be
             stored.  (You must use LAMMPS' "include" command to read this
             file in order to apply these changes.)
          "attribute" specifies the attribute of the atoms that you wish to
             change.  (Eg "type", "charge", "mass". For a list of attributes see
             https://lammps.sandia.gov/doc/set.html)
          a1, a2, ... aM are atom IDs whose attributes will be changed
          i1, i2, ... iM are the monomer indices specifying the relative
             position of each monomer to which atoms a1, a2, ... aM belong
             The i1, i2 arguments are integer indices indicating the (relative)
             location of monomer that each atom belongs to
                0 corresponds to the monomer specified in the "-locations" file
                1 corresponds to the next monomer
                2 corresponds to the following monomer, etc...
	     (For circular polymers, the indices will be wrapped appropriately.)
             These indices are added to the location where you want the
             modification to be made (as specified, for example, in the file
             supplied to the "-locations" argument, or generated automatically).
          A1, A2, ... AM are the new attributes you want these atoms to have
             (For example, the new atom types, new charges, or new masses.)

    -fix-nbody N filename fixname fixID group keyword a1 ... aN i1 ... iN params
             Add fix constraints between specific atoms in the polymer
             using LAMMPS' "fix restrain" or "fix twist" feature.  This allows
             the user a way to apply forces to atoms in the polymer which are
             similar to bond, angle, dihedral, and improper interactions,
             without modifying the bond topology of the system.  Consequently,
             the syntax of this argument is similar to the syntax of the
             -bond, -angle, -dihedral, and -improper arguments.
          The a1, a2, ... aN arguments are strings representing atom names.
          The i1, i2, ... iN arguments are integer offests indicating
             the monomer to which that atom belongs.
          "N" represents the number of atoms effected by the fix
             (2 for bonds, 3 for angles, 4 for dihedrals, impropers, or twists.)
          "filename" is the name of a file where the LAMMPS commands will be
             stored.  (You must use LAMMPS' "include" command to read this file
             in order to enable these new forces during a simulation.)
          "fixname" is the name of the LAMMPS fix you want to use
             (currently only "restrain" and "twist" are supported, but future
              LAMMPS fixes which use similar argument syntax could be used.)
          "fixID" is a string that you will use in your LAMMPS script file
             to identify the fix.  (Typically it is "fxRestrain" or "fxTwist".
             It does not matter as long as it doesn't clash with other fix-IDs.)
          "keyword" is the string specifying what kind of restraint or twist
             motor you want to use.  As of 2020-5-22, for "fix restrain",
             the keywords are "bond", "lbond", "angle", or "dihedral".
             (See https://lammps.sandia.gov/doc/fix_restrain.html)
             For "fix twist", the available keywords are "torque" or "constrain"
             (See https://lammps.sandia.gov/doc/fix_twist.html)
          "group" is the set of atoms you wish to consider when applying this
             fix.  (Atoms outside this set will be ignored.  In LAMMPS, all
             fixes require a group argument.)  Typically it is set to "all".
          "params" is a string containing a list of parameters that describe
             the interaction created by this fix.  For example, if you use
             "fix restrain" to create a "bond", then here you would specify
             the parameters of the spring used to constrain the length of that
             bond (such as the spring stiffness and equilibrium length).
             (In that case, you must surround the two numbers in quotes so
              that the shell does not interpret them as 2 different arguments.)
             If you use "fix twist" to apply torque, then you would specify the
             mangitude of that torque.
          Examples:
           -fix-nbody 2 "fix_link.in" fxBn all restrain bond c c 0 1 "0 10 1.64"
           -fix-nbody 4 "fix_tw.in" fxTw all twist torque b a a b 0 0 1 1 "5.0"

    -circular keyword
          Inform the program that the polymer is circular.
          In order to enable modifications between monomers at opposite
          ends of the (circular) polymer, you must use "-circular yes"
          This allows you to use monomer indices in the polymer which
          may wrap around the polymer.  (Otherwise an error is generated.)
          "keyword" must be one of these choices:
             "no"       The polymer is a linear chain with the two ends
                        not connected.  (default)
             "yes"      The polymer is a circular loop.  (Note: You can also
                        use "connected" instead of "yes".  They are synonymous.)
       

```

## Multiple types of modifications

If you plan to insert multiple *different types* of modifications to the same
polymer, then things are more complicated because you have to figure out
where to place these modifications on the polymer *while preventing overlaps*.
This issue arises, for example, if you have many different kinds of proteins
that you wish to bind to a long DNA polymer (and each of them potentially
occupies a different amount of space on that polymer, or is distributed
in a different way).

In such cases are complicated and you may prefer to write your own script
to generate your own polymer from scratch.  Alternately, you could
create a script which generates the positions of these modifications
and then use "genpoly_modify_lt.py" along with the "-locations" argument
to modify an existing polymer for you.

In this portion of the documentation, I try to explain the procedure that
I use to create custom polymers.
While this procedure is messy, it works for a wide variety of polymers.
No extra scripts are required, but (again) many users my prefer to bypass
this process and create their own polymer builder from scratch.

Recall that all of the modifications made by the "genpoly_modify_lt.py" program
are identical.  So if you have multiple different *types* of modifications
you want to make, you must invoke the "genpoly_modify_lt.py" program
multiple times, once for each type of modification.  Each time, you can
supply it with a different file containing a list of integers
(the -locations file) specifying where you want to modify the polymer.
However you must take care to insure that the different kinds of
modifications are placed at locations on the polymer that do not
overlap with each other.

To help generate the list of non-overlapping integers that you would need, this
program also has some additional features.  For example, you could use the
"-widths", "-locations-random", and "-write-locations" arguments to generate a
list of locations where various types of modifications could be made to a
polymer (without modifying your polymer).
Then you would run the program again to make those modifications.  As mentioned
above, you would do this separately, once for each *type* of modification, each
time reading the corresponding portion of the list of numbers that you already
created earlier.

Alternatively, you could generate the list of integers
in several different stages, using the "-write-occupancy" and "-read-occupancy"
arguments to make sure the resulting integers do not overlap.  This would
allow you to create mixtures of modifications, placing some of them
at regular periodic intervals (or some other distribution you come up with),
and then fitting the remaining modifications around them randomly.

```
    -widths file
             Not all of the modifications need have the same width.  This
             allows you to specify a file containing N positive integers
             (one on each line, where N=number of modifications).  (This is
             typically used when combined with both the -locations-random and
             -write-locations arguments.)

    -write-locations file.txt
             If you are using the -locations-random or -locations-periodic
             arguments, then this argument will cause the program to create a
             file of integers indicating where the modifications will be placed.
             You can use this file to place those objects by invoking this
             program later with the -locations argument.  For example, if the
             first half of the modifications you want to make are of type A
             and the second half of modifications are of type B, you can
             "-locations-random" argument to place these modifications
             randomly on the polymer, and save them to a file using the
             "-write-locations" argument.
             THEN, later you can invoke this program twice to generate the
             LT file that will actually modify the polymer:
             once using the first half of the file to make type A modifications
             once using the scond half of the file to make type B modifications

    -read-occupancy file.txt
             While placing the modifications onto the polymer, you can request
             that some of the sites on the polymer are occupied and unavailable.
             The file argument should be the name of a text file containing a
             list of integers from 0 to N-1, where N is the length of the
             polymer.

    -write-occupancy file.txt
             Write a file containing a list of integers specifying which
             sites on the polymer are occupied (and cannot be used).

    -locations-random-attempts max_attempts
             The algorithm I use with the "-locations-random" argument
             to insert modifications on a polymer is not guaranteed to succeed
             if the user has restricted certain sites using the -read-occupancy
             argument.  It might be simply impossible to fit the desired number
             of modifications along the polymer.  In other cases, it might be
             possible, however the program may fail to do it on the first
             attempt.  This argument allows you to specify how many attempts
             you wish to make before giving up.  A higher number increases the
             chance of success at the cost of slower running times.
```

### Example combining multiple types of modifications on the same polymer.

First we create a polymer.
This example is from the "18bp3p" coarse grained
[DNA example](../examples/coarse_grained/DNA_models/).
```
  genpoly_lt.py \
      -helix 0 \
      -bond Backbone c2 c1 \
      -dihedral Backbone r c2 c2 r 0 0 1 1 \
      -polymer-name 'DNAPolymer' \
      -inherits 'DNAForceField'  \
      -monomer-name 'DNAMonomer' \
      -header 'import "dna_monomer.lt"' \
      < init_crds_polymer_backbone.raw \
      > dna_polymer.lt
```
(Note: The "init_crds_polymer_backbone.raw" file contains 9851 lines containing
 x,y,z coordinates where we want to place the 9851 monomers in the polymer.
 See [genpoly_lt.py documentation](doc_genpoly_lt.md) for details.)

In the next step, we figure out where to distribute 30
modifications (of width 2) at even (ie. periodic) intervals
along the polymer of length 9851.
```
  genpoly_modify_lt.py \
     -length 9851 \
     -locations-periodic 30 0 \
     -width 2 \
     -write-occupancy polymer_occupancy.txt \
     -write-locations locations_twist_motors.txt
```
Now add 514 additional modifications (of width 2) to the polymer, taking
care not to overlap with the modifications you made in the previous step.
```
  genpoly_modify_lt.py \
     -length 9851 \
     -locations-random 514 0 \
     -locations-random-attempts 50 \
     -width 2 \
     -read-occupancy polymer_occupancy.txt \
     -write-occupancy polymer_occupancy.txt \
     -write-locations locations_bends.txt
```
Then we actually modify the polymer.
We place one kind of modification at the positions in the
"locations_twist_motors.txt" file by invoking "genpoly_modify_lt.py" again.
(See the
[twist_motor_supercoiling](examples/coarse_grained/DNA_models/dsDNA_only/twistable_Kratky-Porod/42bp_3particles)
example for details.)
```
  genpoly_modify_lt.py \
     -polymer-name DNAPolymer \
     -length 9851 \
     -locations mod_locations.txt \
     -dihedral Disable r c2 c2 r 0 0 1 1 \
     -fix-nbody 4 "fix_twist.in" fxTw all twist torque r c2 c2 r 0 0 1 1 1.105 \
     -set-atoms 4 "In Types" "type" r c2 c2 r 0 0 1 1 Rm C1m C1m Rm \
    >> dna_polymer.lt
```
*(Note: We used "-width 2" above because this modification changes two different
monomers located at i+0 and i+1, as indicated by the "0 0 1 1" arguments above)*

Finally, we place the other kind of modification at the positions in the
"locations_bends.txt" file by invoking "genpoly_modify_lt.py" again.
```
  genpoly_modify_lt.py \
     -polymer-name DNAPolymer \
     -length 9851 \
     -locations locations_bends.txt \
     -angle Bend c2 c1 c2 0 1 1 \
     -dihedral Disable r c2 c2 r 0 0 1 1 \
     >> dna_polymer.lt
```
At the end of this process, the "dna_polymer.lt" file contains
both types of modifications (twists and bends).


### Warning when using multiple polymers

As of 2020-5-22, for convenience, the "genpoly_lt.py" allows you to create
many different polymers in a single invocation by using the "-cuts" argument.
Each of these polymer objects will be given a different name and length.
In that case, you must examine the LT file and determine the name and length
of each polymer.  Then you can apply this program to each of those polymers,
one at a time.  (In other words, you cannot use a single invocation of this
program to modify the entire system of polymers that you created earlier.
You must invoke it separately for each polymer.)


