interpolate_curve.py
===========

##  Description

"interpolate_curve.py" is a crude program which uses (Catmull-Rom)
cubic spline interpolation to generate a set of evenly spaced
coordinates which lie along smooth a curve specified by the user.

Note: This program is both a stand-alone executable program (that can be run
from the terminal) and a python library.  The former is documented below.
*(The [python API is explained later](#Python-API).)*

## Usage (from the terminal)

```
interpolate_curve.py Ndesired [scale] [alpha] < old_coords.raw > new_coords.raw
```

The old_coords.raw and new_coords.raw are 3-column text files containing
x,y,z coordinates.  The new_coords.raw file in this example will contain
***Ndesired*** coordinates distributed along the path whose control points are
stored in the file "old_coords.raw".  The optional ***scale*** parameter will
cause the resulting coordinates to be multiplied by a constant (default 1).

(The optional ***alpha*** parameter is a number between *0.0* and *1.0*
whose value effects the smoothness of the interpolation.  The default value,
*0.5*, corresponds to a
[centripital Catmull-Rom spline](https://en.wikipedia.org/wiki/Centripetal_Catmull%E2%80%93Rom_spline).)


### Using Moltemplate to trace a polymer along a path

This tool can be useful to help prepare polymer simulations with moltemplate.
There are a variety of algorithms that generate random (space-filling)
curves, and these curves can be used as inputs to the 
[genpoly_lt.py](doc_genpoly_lt.md) program which generates
polymer molecules in moltemplate (.LT) format.
The *genpoly_lt.py* program requires a file containing the
x,y,z coordinates where you want each monomer to go.
An example of its usage can be found
[here.](http://moltemplate.org/images/misc/polymers_follow_a_curve.png)

Unfortunately, most randomly generated curves contain sharp corners,
and the spacing between the control points rarely matches the
distance between monomers in a polymer.
The *interpolate_curve.py* script addresses this issue, allowing
you to resample a curve with the desired scale and resolution.


### Example:

Here we use 
[ndmansfield](https://github.com/jewettaij/ndmansfield)
to generate an initial (space-filling) curve for the shape of the polymer.
Then we use *interpolate_curve.py* to smooth and rescale this curve before
using
[*genpoly_lt.py*](https://github.com/jewettaij/moltemplate/blob/master/doc/doc_genpoly_lt.md)
to build a polymer which can be simulated with the help of
[*moltemplate.sh*](https://github.com/jewettaij/moltemplate)
and
[LAMMPS](http://lammps.sandia.gov).

First step: Create a small jagged version of the curve using *ndmansfield*:
```
ndmansfield -box 53 40 31 -cyclic yes -seed 1  \
            -tsave 200000 -tstop 3200000       \
            > ndmansfield_traj_53x40x31.raw
```
In this example, this program generates a random space-filling curve
(a "Hamiltonian path") which fills a rectangular lattice of size 53 x 40 x 31.
*(This program uses Monte-Carlo to generate a series of increasingly random
polymer shapes as the simulation progresses.)*  The coordinates of the polymer
will be written to a file ("ndmansfield_traj_53x40x31.raw").  This is a
3-column text file (with blank-line delimeters).  It has the following format:

```
x1 y1 z1     #<--1st snapshot in the trajectory
x2 y2 z2
:  :  :
xN yN zN

x1 y1 z1     #<--2nd snapshot in the trajectory
x2 y2 z2
:  :  :
xN yN zN
```
             etc...

*(Notes on "ndmansfield" usage:
We just need to run the lattice simulation long enough to get a random
polymer conformation.  This duration should be long enough.  You can check
long enough by watching the messages from ndmansfield printed to the stderr.
The number of bonds in each direction: x,y,z should be approximately equal.
Incidentally, he coordinates generated ndmansfield are nonnegative integers.)*


We want to extract the last frame of this trajectory.  The length of
the lattice polymer in this example is Nx\*Ny\*Nz = 53\*40\*31 = 65720
but we add 1 because there is a blank line separating each frame.
```
tail -n 65721 ndmansfield_traj_53x40x31.raw > coords.raw
```
It might be convenient to center these coordinates before continuing,
so we do that next:
```
recenter_coords.py 0 0 0 < coords.raw > coords_cen.raw
```
*(Note: The recenter_coords.py script is included with moltemplate.)*

Suppose we want our final polymer to contain 150000 monomers,
run "interpolate_curve.py" this way:
```
interpolate_curve.py 150000 7.99 < coords_cen.raw > coords_smoothed.raw
```
This will generate 150000 points at even intervals along this interpolated
path and rescale them so that the separation distance between monomers
is ~3.50, which equals 7.99\*(65720/150000)\*1.  (Recall that the original
space between the points in the "coords_cen.raw" file is 1.)
Then, use *genpoly_lt.py* to generate a polymer which wraps along the length
of this curve:
```
genpoly_lt.py -monomer-name "Monomer" \
              -polymer-name "Polymer" \
              -bond Backbone C C \
              < coords_smoothed.raw > polymer.lt
```
(Then use moltemplate.sh to convert this file into a
 LAMMPS data file.)



## Python API

```python
def ResampleCurve(x_orig,     # a list or array of points lying along the curve
                  num_points, # number of points you want the new curve to have
                  alpha=0.5)  # optional: the alpha interpolation parameter
```

## Usage example inside python

```python
import numpy as np
N = 10
# Generate a zig-zag curve containing N points
x_orig = np.array([[i, i%2] for i in range(0,N)])

# Now generate a smooth version of this curve:
import moltemplate
x_new = moltemplate.interpolate_curve.ResampleCurve(x_orig, 500, 0.5)
```

Note that there are other free python libraries for curve interpolation, such as
[*scipy.interpolate.interp1d*](https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.interp1d.html).
These libraries are certainly faster (and perhaps more flexible) than this one.
