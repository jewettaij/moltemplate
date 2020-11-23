Long alkane chain example
==============
This example is a simple simulation of a long alkane chain, in a vacuum at room temperature using the COMPASS force field.  The molecule in this example was constructed from monomeric subunits (named "CH2", and "CH3").


#### Images

<img src="images/ch2_ry60.jpg" width=110> <img src="images/plus.svg" height=80> <img src="images/ch3_ry60.jpg" width=110>
<img src="images/rightarrow.svg" height=80> <img src="images/alkane50_t=0_straight.jpg" width=250> <img src="images/rightarrow.svg" height=80> <img src="images/alkane50_t=1ns_equilibrated.jpg" width=150>


### Instructions

1) To build the files which LAMMPS needs, follow the instructions in:
[README_setup.sh](README_setup.sh)

2) To run LAMMPS with these files, follow these instructions:
[README_run.sh](README_run.sh)

*(If you notice a problem with this example, please [report it](../README.md).)*


### Details

The "Alkane50" molecule, as well as the "CH2", and "CH3" monomers it contains, use the COMPASS force-field.  As with all of the COMPASS examples, when we define these molecules, we only specify the atom names, bond list, and coordinates.  We do not have to list the atom charges, angles, dihedrals, or impropers.  The rules for creating atomic charge and angle topology are contained in the ["compass_published.lt"](../../../../moltemplate/force_fields/compass_published.lt) file.  To let moltemplate know that you want to use these rules, define your molecules (and molecular subunits) this way:

```
import "compass_published.lt"
CH2 inherits COMPASS { ... }         # (see "ch2group.lt")
CH3 inherits COMPASS { ... }         # (see "ch3group.lt")
Alkane50 inherits COMPASS { ... }    # (see "alkane50.lt")
```

