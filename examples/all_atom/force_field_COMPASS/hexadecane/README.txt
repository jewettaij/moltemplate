This example is a simple simulation of many long alkane chains (hexadecane) in a
box near the boiling point at atmospheric pressure.  Please read "WARNING.TXT".

NOTE: This particular example uses the COMPASS force-field
      However, moltemplate is not limited to COMPASS.

1) Create the "system.data", "system.in.init", and "system.in.settings"
files which LAMMPS will read by running:

moltemplate.sh system.lt


2) Run LAMMPS in this order:

lmp_mpi -i run.in.npt   # running the simulation at constant pressure
lmp_mpi -i run.in.nvt   # running the simulation at constant temperature

(The name of the LAMMPS executable, eg "lmp_mpi", may vary.)

---- Details ----

The "Hexadecane" molecule, as well as the "CH2", and "CH3" monomers it contains
use the COMPASS force-field.  This means that when we define these molecules,
we only specify the atom names, bond list, and coordinates.
We do not have to list the atom charges, angles, dihedrals, or impropers.
The rules for creating atomic charge and angle topology are contained in
the "compass_published.lt" file created by step 3) above.  The "ch2group.lt",
"ch3group.lt", and "hexadecane.lt" files all refer to "compass_published.lt",
(as well as the "COMPASS" force-field object which it defines).  Excerpt:

import "compass_published.lt"
CH2 inherits COMPASS { ...
CH3 inherits COMPASS { ...
Hexadecane inherits COMPASS { ...

Alternatively, you can manually define a list of angles, dihedrals, and
improper interactions in these files, instead of asking the force-field
to generate them for you.  You can also specify some of the angles and
dihedrals explicitly, and let the force-field handle the rest.
(Many of the examples which come with moltemplate do this.)


