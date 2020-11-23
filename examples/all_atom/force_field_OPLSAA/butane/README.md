Butane example
==============
This example is a simple simulation of many short alkane chains (butane) in a box near the boiling point at atmospheric pressure.  The butane molecule in this example was constructed from monomeric subunits (named "CH2", and "CH3").


#### Images
<img src="images/ch2_ry90.jpg" width=110> <img src="images/plus.svg" height=80> <img src="images/ch3_ry60.jpg" width=110> <img src="images/rightarrow.svg" height=80> <img src="images/butane.jpg" width=150> <img src="images/rightarrow.svg" height=80> <img src="images/initial_configuration_LR.jpg" width=150> <img src="images/rightarrow.svg" height=80> <img src="images/after_pressure_equilibration_LR.jpg" width=150>


### Instructions

1) To build the files which LAMMPS needs, follow the instructions in:
[README_setup.sh](README_setup.sh)

2) To run LAMMPS with these files, follow these instructions:
[README_run.sh](README_run.sh)

*(If you notice a problem with this example, please [report it](../README.md).)*


### Details

The "Butane" molecule, as well as the "CH2", and "CH3" monomers it contains, use the OPLSAA force-field.  As with all of the OPLSAA examples, when we define these molecules, we only specify the atom names, bond list, and coordinates.  We do not have to list the atom charges, angles, dihedrals, or impropers.  The rules for creating atomic charge and angle topology are contained in the ["oplsaa.lt"](../../../../moltemplate/force_fields/oplsaa.lt) file.  To let moltemplate know that you want to use these rules, define your molecules (and molecular subunits) this way:

```
import "oplsaa.lt"
CH2 inherits OPLSAA { ... }      # (see "ch2group.lt")
CH3 inherits OPLSAA { ... }      # (see "ch3group.lt")
Butane inherits OPLSAA { ... }   # (see "butane.lt")
```


### Manual control of bond and angle interactions

It is unlikely that you will need to do this, but if necessary you can customize existing bonds, angles, dihedrals etc. in your molecule (eg. *Butane*), or add new ones (if the force field does not define them).  To do this, edit the corresponding LT file (eg. ["butane.lt"](./moltemplate_files/butane.lt)), and add extra sections to that file (eg. *write("Data Bonds")* or *write("Data Angles")*).  Then add a list of bonded interactions to these sections (containing lines similar to *"\$bond:c7h5 @bond:CustomType \$atom:c7 \$atom:h5"*).  By default, this will override the bond and bonded angular interactions created by the force field.  For more details, read the chapter in the moltemplate manual named "Customizing molecule position and topology".)

