Functionalized Nanotubes
=============================

### Images

<img src="images/graphene_unit_cell.jpg" width=100> <img src="images/plus.svg" height=80> <img src="images/nh2_bbk_occ.jpg" width=60> <img src="images/rightarrow.svg" height=80> <img src="images/graphene_NH2_unit_cell.jpg" width=100> <img src="images/rightarrow.svg" height=80> <img src="images/nanotubes_t=0_bbk.jpg" width=140> <img src="images/rightarrow.svg" height=80> <img src="images/nanotubes_t=100000_bbk.jpg" width=140>


### Details

This example demonstrates a way to build carbon nanotubes with chemical groups attached to the surface at random locations using moltemplate.  In I used the "new" command with [][] brackets to create a 2-D array of graphene unit cells (blue diamonds) which were wrapped around the surface of a cylinder.

A small fraction of unit cells (selected randomly) have a amine group (NH2) attached to one of the carbon atoms (shown above).  This example uses the "new random" command to select randomly from the two different versions of the graphene unit cell (with and without the attached amine group) when filling the 2-D array, as shown in the pictures above.  *(The "new" and "new random" commands are explained in the moltemplate manual.)*


#### There are no carbon-carbon bonds

In the nanotubes, I did not try to connect the carbon atoms together with bonds.  It is possible to build nanotubes with carbon bonds, but this example does not need them.  Instead, the carbon atoms in the nanotubes are rigid (as well as the nitrogen atoms which are directly bonded to them).  However the nanotubes are allowed to move, as are the remaining atoms in the amine groups (the two hydrogen atoms).


#### Force field parameters

The Lennard-Jones parameters for the carbon atoms were taken from this [paper](https://doi.org/10.1016/S0009-2614(01)01127-7).

The 3 atoms in the amine groups use the OPLSAA force field.  The charge of the carbon they are bonded (colored in orange in the pictures above) is modified so that the local structure (4 atoms) is neutral.  In this simulation, these amine groups were assumed to have no other effect on either the shape or the charge of the nearby carbon atoms in the nanotube.


### Requirements

To run this you must have a version of LAMMPS which has been compiled with support for the optional RIGID package. (See the [run.in.nvt](run.in.nvt) file for more details.)  Running at NVT defintely does not require this.

### Notes:

#### Other modeling tools:
If you need explicit bonds between carbon atoms, then you must add them yourself or use a different tool.  Currently (as of 2020-11-30), moltemplate does not generate bonds automatically, although this may change in the future.  The "Nanotube Builder" and "topotools" plugins for for VMD can generate a nanotube with bonds in LAMMPS data format.  You can then convert this data file to .LT format using the ltemplify.py utility and then import it into another .LT file and play with it later.  (In the "cnad-cnt" example, the carbon nanotube was built using "Nanotube Builder" and topotools, and processed with ltemplify.py)

#### Armchair and chiral nanotubes

This is an example of a ["zigzag" nanotube](https://en.wikipedia.org/wiki/Carbon_nanotube#The_zigzag_and_armchair_configurations).  Zigzag nanotubes are easier to build in moltemplate because in a "zigzag" nanotube, the exposed edges at the ends of the nanotube happen to be aligned with the unit cell axes.  If you need to build "armchair" nanotubes or "chiral" *(m,n)* nanotubes, then the ends of the nanotube created using the "new [][]" command will be jagged.  To fix this, you can use moltemplate's "delete" command to remove atoms from the spiky ends of the nanotube.  The "delete" command (and the "new random" command) is discussed in the moltemplate manual.