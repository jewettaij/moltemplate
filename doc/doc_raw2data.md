raw2data.py
===========

## Description

**raw2data.py** replaces the coordinates of a LAMMPS data file with new coordinates.

## Typical usage
```
raw2data.py -atomstyle ATOMSTYLE FILE_OLD.data < COORDS.raw > FILE_NEW.data
```

This will create a new LAMMPS DATA file named "FILE_NEW.data" whose atom
coordinates are copied from the COORDS.raw file, but is otherwise identical
to the original DATA file (eg, "FILE_OLD.data").  The optional
-atomstyle ATOMSTYLE argument tells raw2data.py about the format of the DATA
file.  If not specified, the atom style is "full" by default.


## Arguments

ATOMSTYLE is a quoted string, such as "full" or "hybrid sphere dipole"
indicating the format of the data file.  It can be any of the
[LAMMPS atom styles](https://docs.lammps.org/atom_style.html)
such as "angle", "bond", "charge", "full", "molecular", "dipole", "ellipsoid"
or any hybrid combination of these styles.  (See caveats below.)

FILE_OLD.data
The second argument to raw2data.py is the name of a DATA file you want to read.
raw2data.py will replace the coordinates in the "Atoms" section of this file,
while preserving the rest of the data file.

COORDS.raw is a simple 3-column ASCII file containing the coordinates of the
atoms in your system.  It has a very simple format:
```
-122.28 -19.2293 -7.93705
-121.89 -19.2417 -8.85591
-121.6 -19.2954 -7.20586
-121.59 -20.3273 -2.0079
-122.2 -19.8527 -2.64669
-120.83 -19.7342 -2.2393
  :        :        :
```
The order of the atoms in this file should match the ATOM-ID number in the
first column of the "Atoms" section of the FILE_OLD.data file.
(...I THINK...
 To be on the safe side, use a DATA file with the atoms in sorted order.)

### Exotic atom styles
   When using hybrid atom styles, you must enclose the argument in quotes,
for example: "hybrid sphere dipole"

#### Warning 1

I have not tested using raw2data.py with exotic (non-point-like)
atom styles.  (I suspect that the script will not crash, but dipole
and ellipsoid orientations, and other internal degrees of freedom will
not be updated remain in their initial state.)

Try using [pizza.py](http://pizza.sandia.gov/doc/Manual.html)
instead if you are simulating systems with exotic data types.

#### Warning 2
"raw2data.py" is not a stand-alone script.
Make sure dump2data.py is located in the same directory with raw2data.py.

Note: Although I have not tested it, I suspect many of the other arguments that
work with "dump2data.py", such as "-scale", and "-xyz" also work with raw2data.py.
