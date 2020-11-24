Ice Water (SPC/E)
==============
SPC/E water is arranged in a hexagonal ice (1h) crystal.  The oxygen atoms in 1h ice are arranged in a diamond lattice.  However the orientation of the hydrogen atoms is an additional degree of freedom.  Here I approximated 1h ice as periodic, and provided several different possible choices for the unit cell with different hydrogen bond orientations.

[spce_ice_rect8.lt](moltemplate_files/spce_ice_rect8.lt) (containing 8 water molecules)
[spce_ice_rect16.lt](moltemplate_files/spce_ice_rect16.lt)
[spce_ice_rect32.lt](moltemplate_files/spce_ice_rect32.lt)

The size of the final crystal is a multiple of its unit cell size.

#### Images

Unit Cell (8 water molecules):

<img src="images/ice_rect8_unitcell.png" width=60>

Assembled crystal:

<img src="images/ice_rect8_crystal_3x2x2_LR.jpg" width=150>

The number of water molecules in the simulation and and the simulation box size can be controlled by editing the [system.lt file](moltemplate_files/system.lt).  The simulation contitions can be controlled by editing the [run.in.npt file](run.in.npt).


### Details 

Here I am using the SPC/E water model with long range electrostatics.

To obtain a realistic ice crystal, you may want to run a short simulation to randomize the hydrogen orientations throughout the crystal.  (This can be done at elevated temperatures, imobilizing the oxygen atoms to prevent the ice from melting.  This strategy was not used here.)

### Instructions

More detailed instructions on how to build LAMMPS input files and
run a short simulation are provided in other README files.

1) To build the files which LAMMPS needs, follow the instructions in:
[README_setup.sh](README_setup.sh)

2) To run LAMMPS with these files, follow these instructions:
[README_run.sh](README_run.sh)

*(If you notice a problem with this example, please [report it](../README.md).)*
