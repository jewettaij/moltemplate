MARTINI bilayer preformed (explicit FF parameters)
===============================
This example of the formation of a coarse-grained DPPC lipid-bilayer uses the Martini force-field v2.0 (2013-10), was provided by Saeed Momeni Bashusqeh.  In this version of the example, the lipids are initially arranged in a periodic lattice in the XY plane, pointing in the same direction.  If you prefer, PACKMOL can be used to assemble flat, 2D lipid bilayers with random orientations.  An example using PACKMOL can be found [here](../DPPC_bilayer_formation_PACKMOL/README.md).

### Images
<img src="images/DPPC_martini_LR.jpg" height=110> <img src="images/plus.svg" height=80> <img src="images/water_martini_LR.jpg" width=70> <img src="images/plus.svg" height=80> **PACKMOL** <img src="images/rightarrow.svg" height=80> <img src="images/t=0_bilayer_preformed_GL_LR.jpg" width=170> <img src="images/rightarrow.svg" height=80> <img src="images/t=4ns_bilayer_preformed_GL_LR.jpg" width=170>

The simulation size, number of lipids and water molecules, and initial placement is specified in the [system.lt](./moltemplate_files/system.lt) file.  The simulation contitions can be controlled by editing the [run.in.npt](run.in.npt) and [run.in.nvt](run.in.nvt) files.


### Instructions

More detailed instructions on how to build LAMMPS input files and run a short simulation are provided in other README files.

1) To build the files which LAMMPS needs, follow the instructions in:
[README_setup.sh](README_setup.sh)

2) To run LAMMPS with these files, follow these instructions:
[README_run.sh](README_run.sh)

