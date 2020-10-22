This example is a simple simulation of a long alkane chain,
in a vacuum at room temperature using the OPLSAA force field.

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
lmp_mpi -i run.in.min   # minimize the energy (to avoid atom overlap) before...
lmp_mpi -i run.in.nvt   # running the simulation at constant temperature
```
*(Note: The name of the LAMMPS executable, eg "lmp_mpi", may vary.
See "README_run.sh" for details.)*

### Details

The "Alkane50" molecule, as well as the "CH2", and "CH3" monomers it contains
use the OPLSAA force-field.  This means that when we define these molecules,
we only specify the atom names, bond list, and coordinates.
We do not have to list the atom charges, angles, dihedrals, or impropers.
The rules for creating atomic charge and angle topology are contained in
the "loplsaa.lt" file created by step 3) above.  The "ch2group.lt",
"ch3group.lt", and "alkane50.lt" files all refer to "loplsaa.lt",
(as well as the "OPLSAA" force-field object which it defines).  Excerpt:

```
import "loplsaa.lt"
CH2 inherits OPLSAA { ...
CH3 inherits OPLSAA { ...
Alkane50 inherits OPLSAA { ...
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
(eg. *Alkane50*), or add new ones (if the force field does not define them).
To do this, edit the corresponding LT file (eg. "alkane50.lt"), and add extra
sections to that file (eg *write("Data Bonds")* or *write("Data Angles")*).
Then add a list of bonded interactions to these sections
(containing lines similar to *$bond:cooh1 @bond:004_005 $atom:c7 $atom:o3*).
For more details, read the chapter in the moltemplate manual named
"Customizing molecule position and topology".  It's only a few pages long.)
