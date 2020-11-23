Long alkane chain example
==============
This example is a simple simulation of a long alkane chain, in a vacuum at room temperature using the OPLSAA force field.  The molecule in this example was constructed from monomeric subunits (named "CH2", and "CH3").  In this example, the polymer was initially straight.  However you can use the ["genpoly_lt.py" script](../../../../doc/doc_genpoly_lt.md) to generate .lt files describing long polymers that can wrap around any curve.  (An example of "genpoly_lt.py" usage can be found [here](../../../coarse_grained/DNA_models/dsDNA_only/2strands/3bp_2particles/confined_viral_DNA).)

#### Images

<img src="images/ch2_ry90.jpg" width=110> <img src="images/plus.svg" height=80> <img src="images/ch3_ry60.jpg" width=110>
<img src="images/rightarrow.svg" height=80> <img src="images/alkane50_t=0_straight.jpg" width=250> <img src="images/rightarrow.svg" height=80> <img src="images/alkane50_t=1ns_equilibrated.jpg" width=150>

*Note: This particular example uses the a variant of the OPLSAA force-field suitable for long alkane chains (sometimes called the "LOPLSAA" force field).*


### Instructions

1) To build the files which LAMMPS needs, follow the instructions in:
[README_setup.sh](README_setup.sh)

2) To run LAMMPS with these files, follow these instructions:
[README_run.sh](README_run.sh)

*(If you notice a problem with this example, please [report it](../README.md).)*


### Details

The "Alkane50" molecule, as well as the "CH2", and "CH3" monomers it contains, use the LOPLSAA force-field.  As with all of the OPLSAA examples, when we define these molecules, we only specify the atom names, bond list, and coordinates.  We do not have to list the atom charges, angles, dihedrals, or impropers.  The rules for creating atomic charge and angle topology are contained in the ["loplsaa.lt"](../../../../moltemplate/force_fields/loplsaa.lt) and  ["oplsaa.lt"](../../../../moltemplate/force_fields/oplsaa.lt) files.  To let moltemplate know that you want to use these rules, define your molecules (and molecular subunits) this way:


```
import "loplsaa.lt"
CH2 inherits OPLSAA { ... }         # (see "ch2group.lt")
CH3 inherits OPLSAA { ... }         # (see "ch3group.lt")
Alkane50 inherits OPLSAA { ... }    # (see "alkane50.lt")
```


#### OPLSAA or LOPLSAA"?

There are only a few differences between LOPLSAA and OPLSAA.  The LOPLSAA force field contains a few extra atom types and dihedral interactions which improve the accuracy of long alkane chains.  The ["loplsaa.lt"](../../../../moltemplate/force_fields/loplsaa.lt) file (referenced above) incorporates these extra atom and dihedral types in the existing OPLSAA force field *instead* of creating a new force field named "LOPLSAA".  *(There is no separate "LOPLSAA" force field object.  I apologize if this is confusing.)*  You can mix LOPLSAA and OPLSAA atoms in the same molecule.


### Manual control of bond and angle interactions

If necessary, you can customize existing bonds, angles, dihedrals etc. in your molecule (eg. *Alkane50*), or add new ones (if the force field does not define them).  To do this, edit the corresponding LT file (eg. ["alkane50.lt"](./moltemplate_files/alkane50.lt)), and add extra sections to that file (eg. *write("Data Bonds")* or *write("Data Angles")*).  Then add a list of bonded interactions to these sections (containing lines similar to *"\$bond:c7h5 @bond:CustomType \$atom:c7 \$atom:h5"*).  By default, this will override the bond and bonded angular interactions created by the force field.  For more details, read the chapter in the moltemplate manual named "Customizing molecule position and topology".)

