This example is a simple simulation of many long alkane chains (hexadecane)
in a box near the boiling point at atmospheric pressure.

*Note: This particular example uses the a variant of the OPLSAA force-field
suitable for long alkane chains (sometimes called the "LOPLSAA" force field).
However, moltemplate is not limited to OPLSAA or LOPLSAA.*

### Instructions

1) Create the "system.data", "system.in.init", and "system.in.settings"
files which LAMMPS will read by running:
```
moltemplate.sh system.lt
```
*(See "README_setup.sh" for details.)*

2) Run LAMMPS in this order:
```
lmp_mpi -i run.in.npt   # running the simulation at constant pressure
lmp_mpi -i run.in.nvt   # running the simulation at constant temperature
```
*(Note: The name of the LAMMPS executable, eg "lmp_mpi", may vary.
See "README_run.sh" for details.)*

### Details

The "Hexadecane" molecule, as well as the "CH2", and "CH3" monomers it contains
use the LOPLSAA force-field.  This means that when we define these molecules,
we only specify the atom names, bond list, and coordinates.
We do not have to list the atom charges, angles, dihedrals, or impropers.
The rules for creating atomic charge and angle topology are contained in
the "loplsaa.lt" file created by step 3) above.  The "ch2group.lt",
"ch3group.lt", and "hexadecane.lt" files all refer to "loplsaa.lt",
(as well as the "OPLSAA" force-field objects defined in "loplsaa.lt"
 and "oplsaa.lt").  Excerpt:

```
import "loplsaa.lt"
CH2 inherits OPLSAA { ...
CH3 inherits OPLSAA { ...
Hexadecane inherits OPLSAA { ...
```
#### OPLSAA or LOPLSAA"?

There are only a few differences between LOPLSAA and OPLSAA.
The *import "loplsaa.lt"* line in the exerpt does not create
a new force field named LOPLSAA.  Instead it was implemented in a way
that will modify the existing OPLSAA force field, augmenting it,
adding a few extra atom types and dihedral interactions needed to
improve the accuracy of long alkane chains.
*(I apologize if this is confusing.)*

### Manual control of angle interactions

It is unlikely that you will need to do this, but if necessary
you can customize existing bonds, angles, dihedrals etc. in your molecule
(eg. *Hexadecane*), or add new ones (if the force field does not define them).
To do this, edit the corresponding LT file (eg. "hexadecane.lt"), and add extra
sections to that file (eg *write("Data Bonds")* or *write("Data Angles")*).
Then add a list of bonded interactions to these sections
(containing lines similar to *$bond:cooh1 @bond:004_005 $atom:c7 $atom:o3*).
For more details, read the chapter in the moltemplate manual named
"Customizing molecule position and topology".  It's only a few pages long.)
