dump2data.py
===========

## Description

**dump2data.py** is a tool to extract coordinates from LAMMPS "dump" (trajectory) files.  It was originally designed to convert snapshots from trajectory files into LAMMPS DATA format (for restarting a simulation from where it left off).  However it also reads and writes .XYZ and .RAW (simple 3-column text format) coordinate files.
Although it was written in python, **dump2data.py** is a stand-alone executable intended to be run from the terminal (shell).  *It has not been optimized for speed.*

## Comparison with pizza.py
**dump2data** duplicates some of the tools in
[pizza.py,](http://pizza.sandia.gov/doc/Manual.html).
If you are willing to learn python, **pizza.py** can handle more some dump files which might cause **dump2data.py** to crash.  It includes support for a wider variety of atom styles (eg "atom_style tri").  **dump2data.py** is maintained by the moltemplate developers. **pizza.py** is maintained by the lammps developers.  **pizza.py** may be faster than **dump2data.py** (dump2data.py has not been optimized for speed).


## Arguments

```
   dump2data.py [old_data_file]           \
                [-raw]                    \
                [-xyz]                    \
                [-t time]                 \
                [-tstart ta] [-tstop tb]  \
                [-last]                   \
                [-interval n]             \
                [-type atom_types]        \
                [-id atom_ids]            \
                [-mol mol_ids]            \
                [-multi]                  \
                [-center]                 \
                [-scale x]                \
                [-atomstyle style]        \
                [-xyz-id]                 \
                [-xyz-mol]                \
                [-xyz-type-mol]           \
                < DUMP_FILE > OUTPUT_FILE
```


## Examples

### Example creating RAW (3-column ascii text) coordinate files

   If your LAMMPS dump file is named "traj.lammpstrj", you can
extract the coordinates this way:
```
dump2data.py -raw < traj.lammpstrj > traj.raw
```
The resulting file ("traj.raw") will look like this:
```
-122.28 -19.2293 -7.93705
-121.89 -19.2417 -8.85591
-121.6 -19.2954 -7.20586
  :       :        :

-121.59 -20.3273 -2.0079
-122.2 -19.8527 -2.64669
-120.83 -19.7342 -2.2393
  :       :        :
```
(Note: Blank lines are used to delimit different frames in the trajectory.  If you only want a single frame from the trajectory, you can specify it using the -t or -last arguments.  Alternately you can use the *head* and *tail* unix commands to extract a portion of the trajectory file containing the frame you are interested in beforehand.)

   To limit the output to consider only atoms of a certain type, (for example,
atom types 1,3,5,6, and 7), you can use the *"-type LIST"* argument (for example
"-type 1,3,5-7").  You can also restrict the output to atoms or molecules
with a certain range of ID numbers using the "-id" and "-mol" arguments.
*(These arguments have have no effect if you are creating a ".data" file.
They only work if you are using the '-raw' or '-xyz' arguments.)*


### Example creating XYZ (4-column ascii text) coordinate files

You can extract the coordinates using the .XYZ format this way:
```
dump2data.py -xyz < traj.lammpstrj > traj.xyz
```
This generates a 4-column text file containing the atom-type (first column) followed by the xyz coordinates on each line of each atom (sorted by atomid).  (If you prefer the first column to be something else, you can use the "-xyz-id", "-xyz-mol", and "-xyz-type-mol" arguments instead.)  If there are multiple frames in the trajectory file, it will concatenate them together this way:
```
8192
LAMMPS data from timestep 50000
5 -122.28 -19.2293 -7.93705
3 -121.89 -19.2417 -8.85591
7 -121.6 -19.2954 -7.20586
:   :          :            :
8192
LAMMPS data from timestep 100000
5 -121.59 -20.3273 -2.0079
3 -122.2 -19.8527 -2.64669
7 -120.83 -19.7342 -2.2393
```
(Note: If you want the atom-ID to appear in the first column use "-xyz-id"
       If you want the molecule-ID to appear in the first column use "-xyz-mol"
       If you want the atom-type AND molecule-ID to appear, use "-xyz-type-mol")


## Examples creating DATA files

"dump2data.py" (and "raw2data.py") can also create lammps DATA files.  You must supply them with an existing DATA file containing the correct number of atoms and topology information, and a file containing the coordinates of the atoms.

If your coordinates are stored in a an ordinary 3-column text file ("RAW" file),
you can create the new DATA file this way:
```
raw2data.py -atomstyle ATOM_STYLE data_file < coords.raw  > new_data_file
```
where ATOMSTYLE is a quoted string, such as "full" or "hybrid sphere dipole".
The "-atomstyle ATOM_STYLE" argument is optional.
The default atom_style it is "full".

If your coordinates are stored in a DUMP file (eg "traj.lammpstrj"), 
you can create a new data file this way:
```
dump2data.py -t 10000 data_file < traj.lammpstrj > new_data_file
```
In this example, "10000" is the timestep for the frame you have selected.  You can use *-last* to select the last frame.  If you do not specify the frame you want, multiple data files may be created.  **WARNING: dump2data.py is slow**.  (If you have a long trajectory file, I recommend using the *tail* and *head* unix commands to extract the portion of the trajectory file containing the frame you want before reading it with dump2data.py.  This will be much faster than using the *-t* or *-last* commands.)

(You can use the "-atomstyle" argument with *dump2data.py* as well.)

Creating multiple data files:
The "-multi" command line argument tells "dump2data.py" to generate a new data file for each frame in the trajectory/dump-file.  Those files will have names ending in ".1", ".2", ".3", ...  (If you use the *-interval* argument, frames in the trajectory whose timestep is not a multiple of the interval will be discarded.)  This (probably) occurs automatically whenever the trajectory file contains multiple frames unless you have specified the frame you want (using the *-t* or *-last* arguments)


### Examples using optional command line arguments

If you want to select a particular frame from the trajectory, use:
```
dump2data.py -xyz -t 10000 < traj.lammpstrj > coords.xyz
```
To select the most recent (complete) frame, use:
```
dump2data.py -xyz -last < traj.lammpstrj > coords.xyz
```
(If the last frame is incomplete, this script will attempt to use the previous frame.)

If you want to select multiple frames, but there are too many frames in your trajectory, you can run dump2data.py this way...
```
dump2data.py -xyz -interval 10000 < traj.lammpstrj > traj.xyz
```
...to indicate the desired interval between frames (it must be a multiple of
the save interval).  You can also use "-tstart 500000 and "-tstop 1000000" arguments to limit the output to a particular range of time.  (500000-1000000 in this example).

### Arguments for scaling and centering coordinates

#### -center

This will center the coordinates around the geometric center, so that the average position of the atoms in each frame is located at the origin.  (This script attempts to pay attention to the periodic image flags.  As such, I think this script works with triclinic cells, but I have not tested that feature carefully.)

#### -scale 1.6

This will multiply the coordinates by a constant (eg "1.6")  *(Warning: This argument has not been tested with trajectory files containing periodic image flags: ix iy iz)*

## Limitations

### Speed
The program is slow.  If speed is important to you, you probably should write your own custom script or use pizza.py which might be faster.  (Again, alternatively, you can use the unix *head* and *tail* commands to extract the portion of the trajectory file you are interested in beforehand.)

### Triclinic cells
Support for triclinic cells has been added, but not tested.

### Exotic atom_styles

This script was designed to work with point-like atoms, and it extracts the
x,y,z coordinates (and if present vx,vy,vz velocity)
and it (by default) copies it to the new data being created by this script.

By default, this script assumes you are using "atom_style full".
If you are using some other atom style (eg "hybrid bond dipole"), then you can try to run it this way:
```
dump2data.py -t 10000 \
  -atomstyle "hybrid bond dipole" \
  old_data_file < traj.lammpstrj > new_data_file
```
In general, the -atomstyle argument can be any of the atom styles listed in the
table at:
http://lammps.sandia.gov/doc/atom_style.html
...such as "angle", "bond", "charge", "full", "molecular", "dipole",
"ellipsoid", or any hybrid combination of these styles.
(When using hybrid atom styles, you must enclose the argument in quotes,
for example: "hybrid sphere dipole")

*Warning: I have not tested using dump2data.py with exotic (non-point-like)
atom styles.
I suspect that the script will not crash, but the dipole or ellipsoid
orientations might not be updated and may remain pointing in their
initial directions.
I suspect that "tri", "template", and "body" atom styles will not work at all.*

You can also customize the order columns you want to appear in that file using
-atomstyle ”molid x y z atomid atomtype mux muy muz”.
*(But again, I worry that the mux, muy, muz information in the new data
file might be out of date.)*

Again, try using pizza.py if you are simulating systems with exotic data types.
http://pizza.sandia.gov/doc/Manual.html

I hope this is useful to someone.
